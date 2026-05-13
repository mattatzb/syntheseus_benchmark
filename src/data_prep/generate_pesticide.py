"""Scrape the BCPC Pesticide Compendium for pesticides and their SMILES."""

import logging
import time
from pathlib import Path

import cirpy
import pandas as pd
import pubchempy as pcp
import requests
from bs4 import BeautifulSoup
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

tqdm.pandas()

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

BASE_URL = "http://www.bcpcpesticidecompendium.org/"
INDEX_URL = f"{BASE_URL}index_cn.html"
RAW_CSV = Path("ids.csv")
OUTPUT_CSV = Path("pesticides.csv")
REQUEST_DELAY = 0.5


def _create_session() -> requests.Session:
    """Create an HTTP session with automatic retries on transient errors."""
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def extract_pesticide_names(session: requests.Session) -> list[str]:
    """Extract pesticide page filenames from the A-Z index page.

    Args:
        session: An active requests session.

    Returns:
        List of pesticide page filenames (e.g. ``"abamectin.html"``).
    """
    response = session.get(INDEX_URL, timeout=30)
    response.raise_for_status()

    soup = BeautifulSoup(response.content, "html.parser")
    names: list[str] = []

    for a_tag in soup.find_all("a", href=True):
        names.append(a_tag["href"])

    # Skip the first entry (navigation link, not a pesticide)
    names = names[1:]
    logger.info("Found %d pesticide pages", len(names))
    return names


def _extract_td_text(soup: BeautifulSoup, header_id: str) -> str | None:
    """Return text of the first ``<td>`` matching the given headers attribute."""
    elements = soup.find_all("td", {"headers": header_id})
    if elements:
        return elements[0].text.strip()
    return None


def scrape_pesticide_data(
    session: requests.Session,
    page_names: list[str],
) -> list[dict[str, str | None]]:
    """Scrape CAS, IUPAC name, InChI key, and activity for each pesticide.

    Args:
        session: An active requests session.
        page_names: Filenames returned by :func:`extract_pesticide_names`.

    Returns:
        List of dicts with keys ``name``, ``cas``, ``iupac_name``,
        ``inchi_key``, and ``class``.
    """
    records: list[dict[str, str | None]] = []

    for page_name in tqdm(page_names, desc="Scraping pesticide pages"):
        time.sleep(REQUEST_DELAY)

        try:
            response = session.get(f"{BASE_URL}{page_name}", timeout=30)
            response.raise_for_status()
        except requests.RequestException as exc:
            logger.warning("Failed to fetch %s: %s", page_name, exc)
            continue

        soup = BeautifulSoup(response.text, "html.parser")

        cas = _extract_td_text(soup, "r5")
        activity = _extract_td_text(soup, "r7")
        iupac_name = _extract_td_text(soup, "r3") or "Failed"
        inchi_key = _extract_td_text(soup, "r11") or "Failed"

        # Strip the .html extension to get the compound name
        name = page_name.removesuffix(".html")

        records.append({
            "name": name,
            "cas": cas,
            "iupac_name": iupac_name,
            "inchi_key": inchi_key,
            "class": activity,
        })

    logger.info("Scraped %d records", len(records))
    return records


def get_smiles_from_identifiers(cas: str | None, inchi_key: str | None) -> str:
    """Resolve SMILES from a CAS number (via CIR) or InChI key (via PubChem).

    Args:
        cas: CAS registry number.
        inchi_key: InChI key string.

    Returns:
        SMILES string, or ``"Failed"`` if resolution fails.
    """
    smiles = None

    if cas is not None:
        smiles = cirpy.resolve(cas, "smiles")

    if smiles is None and inchi_key is not None:
        try:
            results = pcp.get_compounds(inchi_key, "inchikey")
            if results:
                smiles = results[0].isomeric_smiles
        except Exception as exc:
            logger.debug("PubChem lookup failed for %s: %s", inchi_key, exc)

    return smiles if smiles is not None else "Failed"


def fetch_or_load_raw(session: requests.Session) -> pd.DataFrame:
    """Load the raw database from disk or scrape it if missing."""
    if RAW_CSV.exists():
        logger.info("Loading cached raw data from %s", RAW_CSV)
        return pd.read_csv(RAW_CSV)

    logger.info("Scraping pesticide index from %s", INDEX_URL)
    page_names = extract_pesticide_names(session)

    logger.info("Scraping individual pesticide pages...")
    records = scrape_pesticide_data(session, page_names)

    df = pd.DataFrame(records)
    df.to_csv(RAW_CSV, index=False)
    logger.info("Saved raw data to %s", RAW_CSV)
    return df


def main() -> None:
    session = _create_session()

    df = fetch_or_load_raw(session)
    df = df.fillna("None")

    logger.info("Resolving SMILES from CAS numbers and InChI keys...")
    df["smiles"] = df.progress_apply(
        lambda row: get_smiles_from_identifiers(row["cas"], row["inchi_key"]),
        axis=1,
    )

    df.to_csv(OUTPUT_CSV, index=False)
    logger.info("Saved %d rows to %s", len(df), OUTPUT_CSV)


if __name__ == "__main__":
    main()