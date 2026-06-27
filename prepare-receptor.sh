SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Select box from pocket-boxes.tsv (column `pocket` is 0-indexed)
if [ -z "$3" ]; then
  echo "Usage: $0 <pdb_file> <output_dir> [box_number]"
  exit 1
fi
POCKET=$3

FOUND=$(awk -F'\t' -v p="$POCKET" 'NR>1 && $1==p { print; exit }' "$SCRIPT_DIR/pocket-boxes.tsv")
if [ -z "$FOUND" ]; then
  echo "Error: no pocket entry with pocket=$POCKET in pocket-boxes.tsv"
  exit 1
fi

center_x=$(echo "$FOUND" | cut -f2)
center_y=$(echo "$FOUND" | cut -f3)
center_z=$(echo "$FOUND" | cut -f4)
side_x=$(echo "$FOUND" | cut -f5)
side_y=$(echo "$FOUND" | cut -f6)
side_z=$(echo "$FOUND" | cut -f7)

mk_prepare_receptor.py --read_pdb "$1" -o "$2" \
  -p -v --box_size "$side_x" "$side_y" "$side_z" --box_center "$center_x" "$center_y" "$center_z"
