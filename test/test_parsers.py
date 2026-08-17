"""Unit tests for the parsing functions.

These cover the pure, table-in/dict-out functions — the ones that read Emu's
output and NCBI's taxonomy. They exist because that is where the real bugs
have been: every case marked REGRESSION below is a bug that actually shipped
and was caught only by running the whole 35-minute pipeline and squinting at
the result. Each now fails in under a second.

Run with:
    python -m pytest test/ -v
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
sys.path.insert(0, str(ROOT / "tools"))

from make_report import (  # noqa: E402
    read_abundance, read_summary, _clean_sample, composition_svg, funnel_svg,
)
from build_emu_db import species_ancestor  # noqa: E402


def write(path, rows, delimiter="\t"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(delimiter.join(map(str, r)) for r in rows) + "\n",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# read_abundance
# ---------------------------------------------------------------------------

def test_lineage_columns_are_not_mistaken_for_samples(tmp_path):
    """REGRESSION: Emu writes the full lineage before the sample columns.

    Treating "every column that isn't tax_id" as a sample turned genus,
    family, order, class, phylum and superkingdom into six phantom barcodes
    that appeared in the report alongside the real ones.
    """
    f = write(tmp_path / "species.tsv", [
        ["tax_id", "species", "genus", "family", "order", "class", "phylum",
         "superkingdom", "barcode01", "barcode02"],
        [654, "Aeromonas veronii", "Aeromonas", "Aeromonadaceae",
         "Aeromonadales", "Gammaproteobacteria", "Pseudomonadota", "Bacteria",
         0.6, 0.4],
        [562, "Escherichia coli", "Escherichia", "Enterobacteriaceae",
         "Enterobacterales", "Gammaproteobacteria", "Pseudomonadota", "Bacteria",
         0.4, 0.6],
    ])
    samples, table = read_abundance(f, "species")
    assert samples == ["barcode01", "barcode02"]
    assert set(table) == {"Aeromonas veronii", "Escherichia coli"}
    assert table["Aeromonas veronii"]["barcode01"] == pytest.approx(0.6)


def test_barcodes_are_ordered_the_same_way_every_run(tmp_path):
    """REGRESSION: Emu's column order is its input directory listing.

    Two runs of the same data produced tables in different orders, and the
    report drew its funnel chart from preprocessing_summary.csv (sorted) and
    its composition charts from this header (not), so the two charts on one
    page listed barcodes differently and could not be read against each other.
    """
    f = write(tmp_path / "s.tsv", [
        ["tax_id", "species", "barcode06", "barcode03", "barcode01"],
        [1, "Escherichia coli", 0.5, 0.2, 0.3],
    ])
    samples, table = read_abundance(f, "species")
    assert samples == ["barcode01", "barcode03", "barcode06"]
    # Re-ordering the names must carry their columns with them.
    assert table["Escherichia coli"]["barcode01"] == pytest.approx(0.3)
    assert table["Escherichia coli"]["barcode03"] == pytest.approx(0.2)
    assert table["Escherichia coli"]["barcode06"] == pytest.approx(0.5)


def test_a_barcode_named_like_a_rank_is_still_excluded_by_name(tmp_path):
    """The lineage filter is by name, so a barcode literally called "class"
    would be dropped. Documented limitation — asserted so a future change to
    value-only detection is a deliberate decision, not an accident."""
    f = write(tmp_path / "s.tsv", [
        ["tax_id", "species", "class", "barcode01"],
        [1, "Escherichia coli", 0.5, 0.5],
    ])
    samples, _ = read_abundance(f, "species")
    assert samples == ["barcode01"]


def test_text_column_with_a_non_rank_name_is_not_a_sample(tmp_path):
    """Name alone is not enough either — a non-rank text column must be
    rejected on its values."""
    f = write(tmp_path / "s.tsv", [
        ["tax_id", "species", "notes", "barcode01"],
        [1, "Escherichia coli", "checked", 1.0],
    ])
    samples, _ = read_abundance(f, "species")
    assert samples == ["barcode01"]


def test_blank_cells_count_as_zero_not_as_errors(tmp_path):
    """Emu leaves a cell empty when a taxon is absent from a barcode."""
    f = write(tmp_path / "s.tsv", [
        ["tax_id", "species", "barcode01", "barcode02"],
        [1, "Escherichia coli", 1.0, ""],
    ])
    _, table = read_abundance(f, "species")
    assert table["Escherichia coli"]["barcode02"] == 0.0


def test_repeated_taxon_rows_are_summed(tmp_path):
    """Combining ranks collapses several species onto one genus."""
    f = write(tmp_path / "g.tsv", [
        ["tax_id", "genus", "barcode01"],
        [1, "Escherichia", 0.3],
        [2, "Escherichia", 0.2],
    ])
    _, table = read_abundance(f, "genus")
    assert table["Escherichia"]["barcode01"] == pytest.approx(0.5)


def test_placeholder_file_yields_nothing(tmp_path):
    """The pipeline writes a '# no non-empty ...' comment when a rank has no
    data. It must parse as empty rather than as a taxon called '#'."""
    f = tmp_path / "s.tsv"
    f.write_text("# no non-empty Emu rel-abundance files available for species\n",
                 encoding="utf-8")
    assert read_abundance(f, "species") == ([], {})


def test_missing_and_empty_files_are_handled(tmp_path):
    assert read_abundance(tmp_path / "nope.tsv", "species") == ([], {})
    empty = tmp_path / "empty.tsv"
    empty.write_text("", encoding="utf-8")
    assert read_abundance(empty, "species") == ([], {})


def test_non_ascii_taxon_names_survive(tmp_path):
    """REGRESSION-adjacent: species names carry accented characters, which is
    why every text read declares utf-8."""
    f = write(tmp_path / "s.tsv", [
        ["tax_id", "species", "barcode01"],
        [1, "Rhodocyclus tenuis Ampuero", 0.5],
        [2, "Zoogloea ramigera Ø-1", 0.5],
    ])
    _, table = read_abundance(f, "species")
    assert "Zoogloea ramigera Ø-1" in table


@pytest.mark.parametrize("raw,expected", [
    ("barcode01_rel-abundance", "barcode01"),
    ("barcode01_filtered", "barcode01"),
    ("barcode01_counts", "barcode01"),
    ("  barcode01  ", "barcode01"),
])
def test_clean_sample_strips_emu_suffixes(raw, expected):
    assert _clean_sample(raw) == expected


# ---------------------------------------------------------------------------
# read_summary
# ---------------------------------------------------------------------------

def test_read_summary_parses_float_formatted_integers(tmp_path):
    """NanoStat writes counts as '993.0'; they must land as ints."""
    f = write(tmp_path / "s.csv", [
        ["barcode", "raw_reads", "filtered_reads"],
        ["barcode01", "993.0", "984.0"],
    ], delimiter=",")
    rows = read_summary(f)
    assert rows[0]["raw"] == 993 and rows[0]["filtered"] == 984


def test_read_summary_tolerates_missing_and_junk_values(tmp_path):
    f = write(tmp_path / "s.csv", [
        ["barcode", "raw_reads", "filtered_reads"],
        ["barcode01", "", "n/a"],
    ], delimiter=",")
    rows = read_summary(f)
    assert rows[0]["raw"] == 0 and rows[0]["filtered"] == 0


# ---------------------------------------------------------------------------
# species_ancestor  (database build)
# ---------------------------------------------------------------------------

# taxid -> (parent, rank).  A small slice shaped like nodes.dmp.
NODES = {
    1: (1, "no rank"),                 # root
    2: (1, "superkingdom"),            # Bacteria
    1224: (2, "phylum"),
    561: (1224, "genus"),              # Escherichia
    562: (561, "species"),             # Escherichia coli
    83334: (562, "no rank"),           # a strain below the species
    99999: (83334, "no rank"),         # deeper still
}


def test_a_species_taxid_maps_to_itself():
    assert species_ancestor(562, NODES) == 562


def test_a_strain_walks_up_to_its_species():
    assert species_ancestor(83334, NODES) == 562
    assert species_ancestor(99999, NODES) == 562


def test_a_genus_has_no_species_ancestor():
    """Emu requires species-level taxids; anything at or above genus is
    dropped rather than silently attached to the wrong rank."""
    assert species_ancestor(561, NODES) is None
    assert species_ancestor(1224, NODES) is None
    assert species_ancestor(2, NODES) is None


def test_unknown_taxid_returns_none():
    assert species_ancestor(4242, NODES) is None


def test_a_cycle_does_not_hang():
    """Defensive: a malformed nodes.dmp must not spin forever."""
    cyclic = {10: (11, "no rank"), 11: (10, "no rank")}
    assert species_ancestor(10, cyclic) is None


# ---------------------------------------------------------------------------
# chart generation — smoke tests, guarding against crashes on edge data
# ---------------------------------------------------------------------------

def test_funnel_handles_a_barcode_with_zero_reads():
    """A barcode that produced nothing must not divide by zero."""
    svg = funnel_svg([{"barcode": "barcode01", "raw": 0, "filtered": 0}])
    assert "<svg" in svg and "barcode01" in svg


def test_composition_handles_no_data():
    svg, legend, ranked = composition_svg([], {}, "species")
    assert "No classification results" in svg
    assert legend == "" and ranked == []


def test_composition_folds_the_tail_into_other():
    """Only eight colour slots exist; a ninth taxon must become "Other"
    rather than getting a generated colour."""
    samples = ["barcode01"]
    table = {f"Species {i}": {"barcode01": 0.1} for i in range(10)}
    svg, legend, _ = composition_svg(samples, table, "species")
    assert "Other" in legend
    assert legend.count("<i class=") == 8  # 7 taxa + Other
