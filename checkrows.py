import sys

def count_onika_rows(file_path):
    """Counts every line in the ONIKA file."""
    try:
        with open(file_path, 'r') as f:
            return sum(1 for line in f)
    except FileNotFoundError:
        print(f"Error: ONIKA file not found at '{file_path}'")
        return 0

def count_dashing_rows(file_path):
    """Counts only data lines (ignoring '#' headers) in the Dashing2 file."""
    try:
        with open(file_path, 'r') as f:
            return sum(1 for line in f if not line.startswith('#'))
    except FileNotFoundError:
        print(f"Error: Dashing2 file not found at '{file_path}'")
        return 0

def main():
    """Main function to execute the script."""
    if len(sys.argv) != 3:
        print("\nUsage: python3 check_rows.py <onika_file_path> <dashing_file_path>")
        sys.exit(1)

    onika_file = sys.argv[1]
    dashing_file = sys.argv[2]
    
    print("--- Verifying Row Counts ---")
    
    onika_rows = count_onika_rows(onika_file)
    dashing_rows = count_dashing_rows(dashing_file)
    
    print(f"ONIKA file ('{onika_file}') has:    {onika_rows} data rows")
    print(f"Dashing2 file ('{dashing_file}') has: {dashing_rows} data rows")
    
    print("-" * 30)
    
    if onika_rows == dashing_rows:
        print("✅ The number of data rows is the same.")
    else:
        print("❌ The number of data rows is DIFFERENT.")
        print("   The main comparison script will fail until these counts match.")

if __name__ == "__main__":
    main()
