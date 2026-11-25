#!/bin/bash

BASE="/mnt/e/AMR Project for Publishing"
LIST_DIR="$BASE/Assembly list"
GENOME_DIR="$BASE/Genomes"

# ============================================================
# Check required programs
# ============================================================

if ! command -v datasets >/dev/null 2>&1; then
    echo "ERROR: NCBI datasets is not installed."
    exit 1
fi

if ! command -v unzip >/dev/null 2>&1; then
    echo "ERROR: unzip is not installed."
    echo "Run: sudo apt update && sudo apt install -y unzip"
    exit 1
fi

# ============================================================
# Create output directories
# ============================================================

mkdir -p "$GENOME_DIR/Animal"
mkdir -p "$GENOME_DIR/Clinical"
mkdir -p "$GENOME_DIR/Environmental"

# ============================================================
# Download function
# ============================================================

download_group() {

    GROUP="$1"
    LIST="$2"
    OUT="$GENOME_DIR/$GROUP"
    LOG="$OUT/download_log.txt"

    mkdir -p "$OUT"

    echo
    echo "============================================================"
    echo "GROUP: $GROUP"
    echo "============================================================"
    echo "Assembly list:"
    echo "$LIST"
    echo
    echo "Output folder:"
    echo "$OUT"
    echo

    if [ ! -f "$LIST" ]; then
        echo "ERROR: Assembly list not found:"
        echo "$LIST"
        return 1
    fi

    total=$(grep -v '^[[:space:]]*$' "$LIST" | wc -l)

    echo "Total assemblies in list: $total"
    echo

    # --------------------------------------------------------
    # Process every assembly
    # --------------------------------------------------------

    while IFS= read -r assembly || [ -n "$assembly" ]; do

        # Remove Windows CR and surrounding whitespace
        assembly=$(echo "$assembly" | tr -d '\r' | xargs)

        # Skip empty lines
        [ -z "$assembly" ] && continue

        # Convert accession to uppercase
        ACCESSION=$(echo "$assembly" | tr '[:lower:]' '[:upper:]')

        echo
        echo "------------------------------------------------------------"
        echo "Processing: $ACCESSION"
        echo "------------------------------------------------------------"

        # ----------------------------------------------------
        # Check whether this assembly already exists
        # ----------------------------------------------------

        if find "$OUT" -maxdepth 1 -type f \
            -iname "${ACCESSION}*.fna" | grep -q .; then

            echo "ALREADY EXISTS: $ACCESSION"
            continue
        fi

        # Temporary ZIP
        ZIP="/tmp/${ACCESSION}.zip"

        # Remove old temporary ZIP if present
        rm -f "$ZIP"

        # ----------------------------------------------------
        # Download genome
        # ----------------------------------------------------

        echo "Downloading $ACCESSION ..."

        if datasets download genome accession "$ACCESSION" \
            --include genome \
            --filename "$ZIP"; then

            echo "Download successful."

        else

            echo "DOWNLOAD FAILED: $ACCESSION"
            echo "DOWNLOAD FAILED: $ACCESSION" >> "$LOG"

            rm -f "$ZIP"

            continue
        fi

        # ----------------------------------------------------
        # Extract all .fna files
        #
        # We intentionally do NOT specify the accession
        # directory. This avoids uppercase/lowercase problems.
        # ----------------------------------------------------

        echo "Extracting .fna file..."

        if unzip -j -o "$ZIP" \
            'ncbi_dataset/data/*/*.fna' \
            -d "$OUT" >/dev/null 2>&1; then

            # ------------------------------------------------
            # Verify extraction
            # ------------------------------------------------

            FNA_COUNT=$(find "$OUT" -maxdepth 1 -type f \
                -iname "${ACCESSION}*.fna" | wc -l)

            if [ "$FNA_COUNT" -gt 0 ]; then

                echo "SUCCESS: $ACCESSION"
                echo "SUCCESS: $ACCESSION" >> "$LOG"

            else

                # Sometimes the filename may not start with
                # exactly the accession, so check any new FNA.
                TOTAL_FNA=$(find "$OUT" -maxdepth 1 \
                    -type f -name "*.fna" | wc -l)

                if [ "$TOTAL_FNA" -gt 0 ]; then

                    echo "FNA extracted successfully: $ACCESSION"
                    echo "SUCCESS: $ACCESSION" >> "$LOG"

                else

                    echo "EXTRACTION FAILED: $ACCESSION"
                    echo "EXTRACTION FAILED: $ACCESSION" >> "$LOG"

                fi
            fi

        else

            echo "EXTRACTION FAILED: $ACCESSION"
            echo "EXTRACTION FAILED: $ACCESSION" >> "$LOG"

        fi

        # ----------------------------------------------------
        # Delete temporary ZIP
        # ----------------------------------------------------

        rm -f "$ZIP"

    done < "$LIST"

    # ========================================================
    # Group summary
    # ========================================================

    downloaded=$(find "$OUT" -maxdepth 1 \
        -type f -name "*.fna" | wc -l)

    echo
    echo "============================================================"
    echo "$GROUP DOWNLOAD COMPLETED"
    echo "============================================================"
    echo "Expected assemblies : $total"
    echo "Downloaded .fna     : $downloaded"
    echo "Output folder       : $OUT"
    echo "============================================================"
    echo
}


# ============================================================
# ANIMAL
# ============================================================

download_group \
    "Animal" \
    "$LIST_DIR/animal_assembly.txt"


# ============================================================
# CLINICAL
# ============================================================

download_group \
    "Clinical" \
    "$LIST_DIR/clinical_assembly.txt"


# ============================================================
# ENVIRONMENTAL
# ============================================================

download_group \
    "Environmental" \
    "$LIST_DIR/environmental_assembly.txt"


# ============================================================
# FINAL SUMMARY
# ============================================================

echo
echo "============================================================"
echo "FINAL GENOME DOWNLOAD SUMMARY"
echo "============================================================"

ANIMAL_COUNT=$(find "$GENOME_DIR/Animal" \
    -maxdepth 1 -type f -name "*.fna" | wc -l)

CLINICAL_COUNT=$(find "$GENOME_DIR/Clinical" \
    -maxdepth 1 -type f -name "*.fna" | wc -l)

ENVIRONMENTAL_COUNT=$(find "$GENOME_DIR/Environmental" \
    -maxdepth 1 -type f -name "*.fna" | wc -l)

echo
echo "Animal genomes        : $ANIMAL_COUNT"
echo "Clinical genomes      : $CLINICAL_COUNT"
echo "Environmental genomes: $ENVIRONMENTAL_COUNT"

echo
echo "Animal folder:"
echo "$GENOME_DIR/Animal"

echo
echo "Clinical folder:"
echo "$GENOME_DIR/Clinical"

echo
echo "Environmental folder:"
echo "$GENOME_DIR/Environmental"

echo
echo "============================================================"
echo "ALL DOWNLOADS COMPLETED"
echo "============================================================"
