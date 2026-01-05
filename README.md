# RNA Secondary Structure Prediction (Nussinov Algorithm)

This project implements the Nussinov algorithm for RNA secondary structure prediction and compares the results with biological reference structures and mfold predictions.

## Project Description

The tool allows users to:
1.  Input an RNA sequence.
2.  Predict its secondary structure using the Nussinov algorithm (maximizing base pairs).
3.  Input a biological reference structure.
4.  Input an mfold predicted structure.
5.  Compare the Nussinov prediction and mfold prediction against the reference using Levenshtein distance.

## Prerequisites

*   Python 3.x
*   `numpy` library

## Installation

1.  Clone the repository or download the files.
2.  Install the required dependencies:

    ```bash
    pip install -r requirements.txt
    ```

## How to Run

Run the main script from the command line:

```bash
python main.py
```

Follow the interactive prompts to enter the RNA sequence and the structures for comparison.

## Project Structure

*   `main.py`: The entry point of the application. Handles user input and displays results.
*   `nussinov_algo.py`: Contains the `NussinovRNA` class which implements the Nussinov algorithm logic.
*   `rna_utils.py`: Contains utility functions, such as `calculate_edit_distance`.
*   `requirements.txt`: List of Python dependencies.
