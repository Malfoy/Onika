import sys
import os

def parse_onika_line(line: str):
    """
    Parses a single line from an ONIKA output file.
    Splits the line by the first tab to separate the description from the data.
    Then, splits the data part by commas.
    
    Args:
        line (str): A single line from the input file.
        
    Returns:
        tuple: A tuple containing:
            - str: The genome description (row header).
            - list[float]: A list of the floating-point values for that row.
    """
    try:
        # Split only on the first tab to correctly handle descriptions with spaces
        description, data_str = line.strip().split('\t', 1)
        
        # Convert the comma-separated string of values into a list of floats
        values = [float(v) for v in data_str.split(',')]
        
        return description, values
    except (ValueError, IndexError) as e:
        # This will catch errors from split not finding a tab, or float() failing on non-numeric data
        raise ValueError(f"Failed to parse line. Ensure it's tab-separated with comma-separated numbers.\nLine: '{line.strip()}'") from e


def analyze_onika_files(file1_path, file2_path):
    """
    Calculates comparison metrics for two ONIKA output files in a memory-efficient way.
    
    Args:
        file1_path (str): Path to the first ONIKA file.
        file2_path (str): Path to the second ONIKA file.
        
    Returns:
        dict: A dictionary containing the calculated difference metrics.
    """
    # Initialize accumulators
    total_absolute_difference = 0.0
    total_squared_difference = 0.0
    element_count = 0
    row_index = 0

    try:
        with open(file1_path, 'r') as f1, open(file2_path, 'r') as f2:
            while True:
                line1 = f1.readline()
                line2 = f2.readline()

                # Stop when both files have been fully read
                if not line1 and not line2:
                    break
                
                # Error if files have a different number of rows
                if not line1 or not line2:
                    raise ValueError("Files have a different number of data rows.")

                row_index += 1
                
                # Parse each line to get description and values
                desc1, values1 = parse_onika_line(line1)
                desc2, values2 = parse_onika_line(line2)

                # --- Validation ---
                if desc1 != desc2:
                    raise ValueError(f"Row description mismatch on line {row_index}.\nFile 1: '{desc1}'\nFile 2: '{desc2}'")
                
                if len(values1) != len(values2):
                    raise ValueError(f"Row {row_index} ('{desc1}') has a different number of values in each file.")

                # --- Calculation ---
                for v1, v2 in zip(values1, values2):
                    raw_diff = v1 - v2
                    total_absolute_difference += abs(raw_diff)
                    total_squared_difference += raw_diff ** 2
                
                element_count += len(values1)

    except FileNotFoundError as e:
        raise FileNotFoundError(f"Error: The file '{e.filename}' was not found.")
    
    if element_count == 0:
        return { "sum_absolute_difference": 0.0, "mean_absolute_difference": 0.0, "sum_squared_difference": 0.0, "mean_squared_difference": 0.0 }
    
    # --- Final Results ---
    results = {
        "sum_absolute_difference": total_absolute_difference,
        "mean_absolute_difference": total_absolute_difference / element_count,
        "sum_squared_difference": total_squared_difference,
        "mean_squared_difference": total_squared_difference / element_count
    }
    return results


def main():
    """
    Main function to execute the script.
    """
    if len(sys.argv) != 3:
        print("Usage: python3 your_script_name.py <file1_path> <file2_path>")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    try:
        metrics = analyze_onika_files(file1_path, file2_path)
        
        print("\n--- Absolute Difference ---")
        print(f"Sum of Absolute Differences:      {metrics['sum_absolute_difference']}")
        print(f"Mean of Absolute Differences:     {metrics['mean_absolute_difference']}")
        
        print("\n--- Squared Difference ---")
        print(f"Sum of Squared Differences:       {metrics['sum_squared_difference']}")
        print(f"Mean of Squared Differences:      {metrics['mean_squared_difference']}")

    except (FileNotFoundError, ValueError) as e:
        print(f"\nAn error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()