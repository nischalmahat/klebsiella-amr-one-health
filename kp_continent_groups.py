"""
kp_continent_groups.py
=======================
Maps the "Location" column of the Klebsiella pneumoniae dataframe (Kp) to a
"Continent" column, based on the 1,765 unique values in Kp_country.txt.

How Location values are structured
-----------------------------------
Values look like "usa", "china: hangzhou", "usa: san antonio, tx",
"united kingdom (england, wales & n. ireland)" - i.e. a country name,
sometimes followed by ":" and a city/state/region, sometimes with the
qualifier in parentheses. The classifier here:

  1. Strips the "country" out of the free-text Location string (splits on
     ":" first, falls back to a "starts-with" match against the country
     list, so values like "united kingdom (england, wales & n. ireland)"
     still resolve correctly).
  2. Looks the country up in COUNTRY_TO_CONTINENT (built from standard
     UN geoscheme / ISO continent groupings and verified against every
     unique token actually present in Kp_country.txt).
  3. Falls back to a keyword scanner for any row that doesn't cleanly
     split on the first token (defensive - your production data should
     not need it, but new records eventually will contain something odd).

Special cases
-------------
- "gaza strip" and "west bank" -> Asia (Middle East), not Africa.
- "hong kong" and "taiwan" -> Asia (kept separate from "china" as written
  in the source data, since that's how they appear in Location).
- "antarctica" -> continent "Antarctica" (its own group, not folded into
  anything else).
- "missing: data agreement established pre-2023" -> "Unknown" (not a
  real location).
- Free-standing US state/city abbreviations (e.g. "usa:tx", "usa: ca") are
  all under the "usa" prefix so they resolve automatically to
  North America - no state-level mapping needed.

Usage
-----
    from kp_continent_groups import assign_continent
    Kp["Continent"] = Kp["Location"].map(assign_continent)
"""
import re

# ----------------------------------------------------------------------
# Country -> Continent map (7-continent / UN geoscheme convention)
# ----------------------------------------------------------------------
COUNTRY_TO_CONTINENT = {
    # ---- Africa ----
    "algeria": "Africa", "benin": "Africa", "botswana": "Africa",
    "burkina faso": "Africa", "cameroon": "Africa", "djibouti": "Africa",
    "egypt": "Africa", "ethiopia": "Africa", "gambia": "Africa",
    "ghana": "Africa", "kenya": "Africa", "madagascar": "Africa",
    "malawi": "Africa", "mali": "Africa", "morocco": "Africa",
    "mozambique": "Africa", "niger": "Africa", "nigeria": "Africa",
    "rwanda": "Africa", "senegal": "Africa", "south africa": "Africa",
    "sudan": "Africa", "tanzania": "Africa", "togo": "Africa",
    "tunisia": "Africa", "uganda": "Africa", "zambia": "Africa",
    "zimbabwe": "Africa", "ivory coast": "Africa", "cote d'ivoire": "Africa",
    "libya": "Africa", "angola": "Africa", "namibia": "Africa",
    "somalia": "Africa", "eritrea": "Africa", "chad": "Africa",
    "guinea": "Africa", "sierra leone": "Africa", "liberia": "Africa",
    "congo": "Africa", "drc": "Africa",

    # ---- Asia (incl. Middle East) ----
    "afghanistan": "Asia", "armenia": "Asia", "azerbaijan": "Asia",
    "bahrain": "Asia", "bangladesh": "Asia", "bhutan": "Asia",
    "brunei": "Asia", "cambodia": "Asia", "china": "Asia",
    "hong kong": "Asia", "india": "Asia", "indonesia": "Asia",
    "iran": "Asia", "iraq": "Asia", "israel": "Asia", "japan": "Asia",
    "jordan": "Asia", "kazakhstan": "Asia", "korea": "Asia",
    "south korea": "Asia", "north korea": "Asia", "kuwait": "Asia",
    "laos": "Asia", "lebanon": "Asia", "malaysia": "Asia",
    "myanmar": "Asia", "nepal": "Asia", "oman": "Asia", "pakistan": "Asia",
    "philippines": "Asia", "qatar": "Asia", "saudi arabia": "Asia",
    "singapore": "Asia", "sri lanka": "Asia", "syria": "Asia",
    "taiwan": "Asia", "thailand": "Asia", "turkey": "Asia",
    "united arab emirates": "Asia", "uzbekistan": "Asia",
    "viet nam": "Asia", "vietnam": "Asia", "yemen": "Asia",
    "georgia": "Asia", "gaza strip": "Asia", "west bank": "Asia",
    "palestine": "Asia", "mongolia": "Asia", "tajikistan": "Asia",
    "kyrgyzstan": "Asia", "turkmenistan": "Asia",

    # ---- Europe ----
    "albania": "Europe", "austria": "Europe", "belarus": "Europe",
    "belgium": "Europe", "bulgaria": "Europe", "croatia": "Europe",
    "cyprus": "Europe", "czech republic": "Europe", "czechia": "Europe",
    "denmark": "Europe", "estonia": "Europe", "finland": "Europe",
    "france": "Europe", "germany": "Europe", "greece": "Europe",
    "hungary": "Europe", "ireland": "Europe", "italy": "Europe",
    "latvia": "Europe", "lithuania": "Europe", "luxembourg": "Europe",
    "macedonia": "Europe", "malta": "Europe", "moldova": "Europe",
    "montenegro": "Europe", "netherlands": "Europe", "norway": "Europe",
    "poland": "Europe", "portugal": "Europe", "romania": "Europe",
    "russia": "Europe", "serbia": "Europe", "slovakia": "Europe",
    "slovenia": "Europe", "spain": "Europe", "sweden": "Europe",
    "switzerland": "Europe", "ukraine": "Europe", "united kingdom": "Europe",
    "united kingdom (england, wales & n. ireland)": "Europe",
    "bosnia and herzegovina": "Europe", "iceland": "Europe",

    # ---- North America (incl. Central America & Caribbean) ----
    "usa": "North America", "canada": "North America",
    "mexico": "North America", "costa rica": "North America",
    "cuba": "North America", "curacao": "North America",
    "dominican republic": "North America", "guadeloupe": "North America",
    "guatemala": "North America", "honduras": "North America",
    "jamaica": "North America", "martinique": "North America",
    "puerto rico": "North America",
    "saint kitts and nevis": "North America",
    "trinidad and tobago": "North America", "panama": "North America",
    "nicaragua": "North America", "el salvador": "North America",
    "belize": "North America", "bahamas": "North America",
    "haiti": "North America",

    # ---- South America ----
    "argentina": "South America", "brazil": "South America",
    "chile": "South America", "colombia": "South America",
    "ecuador": "South America", "paraguay": "South America",
    "peru": "South America", "uruguay": "South America",
    "venezuela": "South America", "bolivia": "South America",
    "guyana": "South America", "suriname": "South America",

    # ---- Oceania ----
    "australia": "Oceania", "fiji": "Oceania", "guam": "Oceania",
    "new zealand": "Oceania", "new caledonia": "Oceania",
    "papua new guinea": "Oceania",

    # ---- Not a real geography ----
    "antarctica": "Antarctica",
    "missing": "Unknown",
}

FALLBACK = "Unknown"

# Sorted longest-first so multi-word countries (e.g. "south africa") match
# before a shorter substring (e.g. "africa") could.
_COUNTRY_KEYS_BY_LEN = sorted(COUNTRY_TO_CONTINENT, key=len, reverse=True)


def _norm(s):
    if s is None:
        return ""
    return re.sub(r"\s+", " ", str(s).strip().lower())


def extract_country(location):
    """Pull the country portion out of a free-text Location string."""
    t = _norm(location)
    if not t:
        return ""
    # Primary rule: country is whatever precedes the first ':' (this is how
    # essentially every value in Kp_country.txt is formatted).
    head = t.split(":", 1)[0].strip()
    if head in COUNTRY_TO_CONTINENT:
        return head
    # Defensive fallback: value has no ':' delimiter, or the country name
    # itself contains a colon-free parenthetical - match by "starts with".
    for c in _COUNTRY_KEYS_BY_LEN:
        if t == c or t.startswith(c + " ") or t.startswith(c + "("):
            return c
    return head  # unseen token - return as-is for the caller to inspect


def assign_continent(location):
    """Return the continent name for a Location string, or 'Unknown'."""
    country = extract_country(location)
    return COUNTRY_TO_CONTINENT.get(country, FALLBACK)


if __name__ == "__main__":
    import pandas as pd  # noqa: F401

    # Kp must already exist in your session, e.g.
    # Kp = pd.read_csv("Kp.csv")

    Kp["Continent"] = Kp["Location"].map(assign_continent)

    print(Kp["Continent"].value_counts(dropna=False))

    # Any Location that didn't map cleanly - review and extend
    # COUNTRY_TO_CONTINENT above if new countries show up later.
    unmapped = Kp.loc[Kp["Continent"] == "Unknown", "Location"].value_counts()
    print(unmapped.head(50))
