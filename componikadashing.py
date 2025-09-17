import sys

def analyze_files_simplified(onika_path, dashing_path):
    """
    Compares ONIKA and Dashing2 files based on simplified rules.
    
    Args:
        onika_path (str): Path to the ONIKA matrix file.
        dashing_path (str): Path to the Dashing2 matrix file.
        
    Returns:
        dict: A dictionary containing the calculated metrics.
    """
    try:
        with open(onika_path, 'r') as f_onika, open(dashing_path, 'r') as f_dashing:
            # Rule 1: Ignore the first three lines of the Dashing2 output
            for _ in range(3):
                f_dashing.readline()

            # Initialize accumulators
            total_absolute_difference = 0.0
            total_squared_difference = 0.0
            row_index = 0

            # Process both files line-by-line
            while True:
                line_onika = f_onika.readline()
                line_dashing = f_dashing.readline()

                # Stop when both files end
                if not line_onika and not line_dashing:
                    break
                
                # Error if one file ends before the other
                if not line_onika or not line_dashing:
                    raise ValueError("Files have a different number of data rows after skipping the Dashing2 header.")

                row_index += 1
                
                # Rule 2: Ignore the first element of each line
                try:
                    # Parse ONIKA: split by tab, take data, split by comma
                    onika_data_str = line_onika.strip().split('\t', 1)[1]
                    onika_values = [float(v) for v in onika_data_str.split(',')]
                    
                    # Parse Dashing2: split by tab, take all but the first element
                    dashing_parts = line_dashing.strip().split('\t')
                    dashing_values = [float(v) for v in dashing_parts[1:]]
                except (IndexError, ValueError):
                    raise ValueError(f"Could not parse data on line {row_index}. Please check file formatting.")

                # Validate that the remaining data has the same number of columns
                if len(onika_values) != len(dashing_values):
                    raise ValueError(f"Row {row_index} has a different number of values after ignoring the first column.")

                # Calculate differences for the corresponding rows
                for v1, v2 in zip(onika_values, dashing_values):
                    raw_diff = v1 - v2
                    total_absolute_difference += abs(raw_diff)
                    total_squared_difference += raw_diff ** 2
            
            if row_index == 0:
                print("Warning: No data rows found to compare.")

            return {
                "sum_absolute_difference": total_absolute_difference,
                "sum_squared_difference": total_squared_difference
            }

    except FileNotFoundError as e:
        raise FileNotFoundError(f"Error: A file was not found. Details: {e}")


def main():
    """Main function to execute the script."""
    if len(sys.argv) != 3:
        print("\nUsage: python3 simplified_compare.py <onika_file_path> <dashing_file_path>")
        sys.exit(1)

    onika_file = sys.argv[1]
    dashing_file = sys.argv[2]

    try:
        metrics = analyze_files_simplified(onika_file, dashing_file)
        
        print("\n--- Comparison Results ---")
        print(f"Sum of Absolute Difference: {metrics['sum_absolute_difference']}")
        print(f"Sum of Squared Difference:  {metrics['sum_squared_difference']}")

    except (FileNotFoundError, ValueError) as e:
        print(f"\nAn error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()