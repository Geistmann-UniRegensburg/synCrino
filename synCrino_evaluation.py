import pandas as pd
# Suppress pandas warning
pd.set_option('future.no_silent_downcasting', True)
import numpy as np

import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import linear_sum_assignment

from tqdm import tqdm

def load_database(file_path_database):
    """Load and clean the reference database from an standardized database Excel file."""
    database = pd.read_excel(file_path_database)

    # Extract numeric shift data, clean whitespace, and sort descending per column
    df_numeric = database.iloc[3:, 1:].copy()
    
    for col in df_numeric.columns:
        cleaned_col = pd.to_numeric(df_numeric[col].astype(str).str.replace(r"\s+", "", regex=True), errors='coerce')
        sorted_values = cleaned_col.sort_values(ascending=False, na_position='last').values
        df_numeric[col] = sorted_values

    # Re-insert the cleaned and sorted numerical values back into the original dataframe structure    
    database.iloc[3:, 1:] = df_numeric
    
    return database

def load_sample(file_path_sample):
    """Load and parse experimental 13C NMR shifts from a standardized sample Excel file."""
    sample = pd.read_excel(file_path_sample)
    # Convert to numeric, round to 1 decimal place, drop NaNs, and sort descending
    sample_sorted = pd.to_numeric(sample["13 C [ppm]"], errors='coerce').round(1).dropna().sort_values(ascending=False).values
    return sample_sorted

def best_shift_differences_brute_force(sample, db_values):
    """Brute-force baseline implementation for validation: uses a dense, single-pass grid search (1001 points) combined with the Hungarian algorithm 
    to find the optimal shift correction and RMSD, prioritizing exhaustive precision over computational speed."""
    sample = np.asarray(sample, dtype=float)
    db_values = np.asarray(db_values, dtype=float)

    if len(sample)==len(db_values):
        # For equal-length sets (same carbon count), calculate the uniform systematic 
        # shift correction via the average difference
        best_shift=np.average(sample-db_values)
        shifted_sample = sample - best_shift

        # Compute RMSD (Root Mean Square Deviation) as the primary ranking metric
        differences = (shifted_sample-db_values)**2
        rmsd = np.sqrt(np.mean(differences))
        unmatched_carbons = np.array([])
        
    else:
        # Assign 'a' to the smaller array and 'b' to the larger array for bipartite matching dimensions
        if(len(sample)) < len(db_values):
            a = sample
            b = db_values
            sample_is_smaller = True
        else:
            a = db_values
            b = sample
            sample_is_smaller = False
        
                
        # High-resolution brute-force grid search: evaluate 1001 evenly spaced candidate 
        # shifts across the full interval [-5, 5] ppm to establish an exhaustive baseline
        test_shifts = np.linspace(-5, 5, 1001)
        best_rmsd = float('inf')
        best_shift = 0.0
    
        for s in test_shifts:
            diff_matrix = ((a[:, np.newaxis] - s) - b)**2
            row_ind, col_ind = linear_sum_assignment(diff_matrix)
            score = np.sqrt(np.mean((a[row_ind] - s - b[col_ind])**2))
    
            if score < best_rmsd:
                best_rmsd = score
                best_shift = s
            
        # Final assignment using the brute-force optimized shift
        final_diff_matrix = ((a[:, np.newaxis] - best_shift) - b) ** 2
        a_indices, b_indices = linear_sum_assignment(final_diff_matrix)

        rmsd = np.sqrt(np.mean((a[a_indices] - best_shift - b[b_indices])**2))
  
        # Identify unmatched signals in the larger dataset
        matched_carbons_indices = set(b_indices)
        all_b_indices = set(range(len(b)))
        unmatched_carbons_indices = list(all_b_indices - matched_carbons_indices)
        unmatched_carbons = b[unmatched_carbons_indices]
       
        # Invert the shift sign if sample and database roles were swapped during optimization
        if not sample_is_smaller:
            best_shift = -best_shift
        
    return rmsd,best_shift, unmatched_carbons

def comparison_similar_compound_bruteforce(database, sample, exact_carbon_count = True):
    """Brute-force validation wrapper: iterates through the entire reference database,
    applies the  brute-force matching method for each compound, returns match metrics."""
    # Prepare database: extract numerical shift values starting from row 3, column 1
    database_sorted = database.copy()
    database_sorted = database_sorted.iloc[3:,1:]
    
    sample_sorted = np.asarray(sample, dtype=np.float64)

    if len(sample_sorted) == 0:
        raise ValueError("sample has no peaks.")
    
    n_sample_carbons = len(sample_sorted)

    results = []

    # Iterate through each reference compound column in the database with a progress bar
    for col in tqdm(database_sorted.columns, desc="Comparing compounds"):
        db_values = pd.to_numeric(database_sorted[col], errors='coerce').dropna().values
        if len(db_values) == 0:
            continue

        if exact_carbon_count:
            if len(db_values) != n_sample_carbons:
                continue

        # Compute optimal shift correction, RMSD, and unmatched signals using brute force
        rmsd, best_shift, unmatched_carbons = best_shift_differences_brute_force(sample_sorted, db_values)
        
        # Store match metrics and metadata for validation and benchmarking
        results.append({
            'compound': col,
            'RMSD': round(rmsd, 4),
            'shift': round(best_shift, 2),
            'solvent': database.loc[0, col],
            'carbons': len(db_values),
            'unmatched carbons': np.sort(unmatched_carbons)[::-1],
            'Author + Year': database.loc[1, col],
            'DOI': database.loc[2, col]
        })

    return results

def best_shift_differences_adaptive_grid_search(sample, db_values):
    """Adaptive grid-search baseline: utilizes a 3-step coarse-to-fine resolution strategy combined 
    with the Hungarian algorithm to determine the optimal shift correction, RMSD, and unmatched signals."""
    sample = np.asarray(sample, dtype=float)
    db_values = np.asarray(db_values, dtype=float)

    if len(sample) == len(db_values):
        # For equal-length sets (same carbon count), calculate the uniform systematic 
        # shift correction via the average difference
        best_shift = np.average(sample - db_values)
        differences = (sample - best_shift - db_values) ** 2
        rmsd = np.sqrt(np.mean(differences))
        unmatched_carbons = np.array([])
        return rmsd, best_shift, unmatched_carbons
        
    else:
        # Assign 'a' to the smaller array and 'b' to the larger array for bipartite matching dimensions
        if len(sample) < len(db_values):
            a = sample
            b = db_values
            sample_is_smaller = True
        else:
            a = db_values
            b = sample
            sample_is_smaller = False
        
        best_rmsd = float('inf')
        best_shift = 0.0

        # To handle datasets of unequal length, we treat the signal alignment as an 
        # optimal bipartite matching problem using an adaptive 3-step grid search:

        # 1. Coarse search: Test a wide range of potential systematic shifts (-5.0 to +5.0 ppm)
        shifts_coarse = np.linspace(-5, 5, 51) 
        for s in shifts_coarse:
            diff_matrix = ((a[:, np.newaxis] - s) - b) ** 2
            row_ind, col_ind = linear_sum_assignment(diff_matrix)
            score = np.sqrt(np.mean((a[row_ind] - s - b[col_ind]) ** 2))

            if score < best_rmsd:
                best_rmsd = score
                best_shift = s

        # 2. Medium resolution search: Refine search window around the best coarse shift
        shifts_medium = np.linspace(best_shift - 0.4, best_shift + 0.4, 17)
        for s in shifts_medium:
            # Ensure the evaluated shift stays within valid boundaries
            if s < -5 or s > 5: 
                continue 
        
            diff_matrix = ((a[:, np.newaxis] - s) - b) ** 2
            row_ind, col_ind = linear_sum_assignment(diff_matrix)
            score = np.sqrt(np.mean((a[row_ind] - s - b[col_ind]) ** 2))

            if score < best_rmsd:
                best_rmsd = score
                best_shift = s

        # 3. Fine resolution search: Final high-precision optimization around the best medium shift
        shifts_fine = np.linspace(best_shift - 0.1, best_shift + 0.1, 21)
        for s in shifts_fine:
            if s < -5 or s > 5: 
                continue

            diff_matrix = ((a[:, np.newaxis] - s) - b) ** 2
            row_ind, col_ind = linear_sum_assignment(diff_matrix)
            score = np.sqrt(np.mean((a[row_ind] - s - b[col_ind]) ** 2))

            if score < best_rmsd:
                best_rmsd = score
                best_shift = s
            
        # Final assignment using the optimized shift
        final_diff_matrix = ((a[:, np.newaxis] - best_shift) - b) ** 2
        a_indices, b_indices = linear_sum_assignment(final_diff_matrix)

        rmsd = np.sqrt(np.mean((a[a_indices] - best_shift - b[b_indices]) ** 2))

        # Identify unmatched signals in the larger dataset
        matched_carbons_indices = set(b_indices)
        all_b_indices = set(range(len(b)))
        unmatched_carbons_indices = list(all_b_indices - matched_carbons_indices)
        unmatched_carbons = b[unmatched_carbons_indices]
        
        # Invert the shift sign if sample and database roles were swapped during optimization
        if not sample_is_smaller:
            best_shift = -best_shift
        
        return rmsd, best_shift, unmatched_carbons

def comparison_similar_compound_ags(database, sample, exact_carbon_count = True):
    """Adaptive grid search validation wrapper: iterates through the entire reference database,
    applies the adaptive 3-step grid-search matching method for each compound, and returns match metrics."""
    # Prepare database: extract numerical shift values starting from row 3, column 1
    database_sorted = database.copy()
    database_sorted = database_sorted.iloc[3:,1:]
    
    sample_sorted = np.asarray(sample, dtype=np.float64)

    if len(sample_sorted) == 0:
        raise ValueError("sample has no peaks.")
    
    n_sample_carbons = len(sample_sorted)

    results = []

    # Iterate through each reference compound column in the database with a progress bar
    for col in tqdm(database_sorted.columns, desc="Comparing compounds"):
        db_values = pd.to_numeric(database_sorted[col], errors='coerce').dropna().values
        if len(db_values) == 0:
            continue

        if exact_carbon_count:
            if len(db_values) != n_sample_carbons:
                continue
    
        # Compute optimal shift correction, RMSD, and unmatched signals using the adaptive grid search method
        rmsd, best_shift, unmatched_carbons = best_shift_differences_adaptive_grid_search(sample_sorted, db_values)
        
        # Store match metrics and metadata for validation and benchmarking
        results.append({
            'compound': col,
            'RMSD': round(rmsd, 4),
            'shift': round(best_shift, 2),
            'solvent': database.loc[0, col],
            'carbons': len(db_values),
            'unmatched carbons': np.sort(unmatched_carbons)[::-1],
            'Author + Year': database.loc[1, col],
            'DOI': database.loc[2, col]
        })

    return results

def dp_monotonic_match(a, b, shift):
    """Dynamic programming helper for monotonic signal alignment: finds the optimal subset-matching of array 'a' 
    into array 'b' under a fixed shift correction, preserving the strict sequential (monotonic) order of chemical shifts."""

    n, m = len(a), len(b)
    INF = np.inf
    # DP matrix size is (n+1) x (m+1) to accommodate the 0-element boundary condition
    dp = np.full((n + 1, m + 1), INF)
    dp[0, :] = 0.0
    back = np.zeros((n + 1, m + 1), dtype=bool)  # Tracks traceback path (True = matched here)

    for i in range(1, n + 1):
        for j in range(i, m + 1):
            # Option 1: Match element a[i-1] with b[j-1] under the given shift
            cost_match = dp[i - 1, j - 1] + (a[i - 1] - shift - b[j - 1]) ** 2
            # Option 2: Skip element b[j-1] if maintaining monotonicity allows it
            cost_skip = dp[i, j - 1] if j > i else INF
            if cost_match <= cost_skip:
                dp[i, j] = cost_match
                back[i, j] = True
            else:
                dp[i, j] = cost_skip
                back[i, j] = False

    # Backtrack through the DP matrix to retrieve the optimal matched indices in b
    i, j = n, m
    matched_b = []
    while i > 0:
        if back[i, j]:
            matched_b.append(j - 1)
            i -= 1
            j -= 1
        else:
            j -= 1
    matched_b.reverse()

    return dp[n, m], matched_b


def best_shift_dp(sample, db_values, shift_bounds=(-5, 5), grid_points=51):
    """Dynamic programming baseline for validation: combines a dense shift grid search with monotonic dynamic programming
    to evaluate shift-corrected alignments while preserving spectral order."""

    sample = np.sort(np.asarray(sample, dtype=float))
    db_values = np.sort(np.asarray(db_values, dtype=float))

    if len(sample) == len(db_values):
        # For equal-length sets (same carbon count), calculate the uniform systematic 
        # shift correction via the average difference
        best_shift = np.mean(sample - db_values)
        rmsd = np.sqrt(np.mean((sample - best_shift - db_values) ** 2))
        return rmsd, best_shift, np.array([])

    # Assign 'a' to the smaller array and 'b' to the larger array for alignment dimensions
    if len(sample) < len(db_values):
        a, b = sample, db_values
        sample_is_smaller = True
    else:
        a, b = db_values, sample
        sample_is_smaller = False

    lo, hi = shift_bounds
    grid = np.linspace(lo, hi, grid_points)

    best_rmsd = np.inf
    best_shift = None
    best_matched_b = None

    # Iterate over candidate shifts using the dynamic programming alignment core
    for s in grid:
        cost, matched_b = dp_monotonic_match(a, b, s)
        
        # Refine the shift using the exact closed-form mean difference of the matched subset
        b_sub = b[matched_b]
        refined_shift = np.clip(np.mean(a - b_sub), lo, hi)
        rmsd = np.sqrt(np.mean((a - refined_shift - b_sub) ** 2))

        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_shift = refined_shift
            best_matched_b = matched_b

    # Identify unmatched signals in the larger dataset
    all_b_idx = set(range(len(b)))
    unmatched_idx = list(all_b_idx - set(best_matched_b))
    unmatched_carbons = b[unmatched_idx]

    # Invert the shift sign if sample and database roles were swapped during optimization
    final_shift = best_shift if sample_is_smaller else -best_shift
    return best_rmsd, final_shift, unmatched_carbons

def compute_all_matches(sample, database,
                         shift_bounds=(-5, 5), grid_points=51, exact_carbon_count = True):
    """Database-wide matching wrapper for the dynamic programming (monotonic) approach: 
    iterates through all reference compounds using a progress bar, computes optimal 
    alignments, and compiles comprehensive metrics into a list of dictionaries."""
    sample = np.asarray(sample, dtype=float)
    results = []
    n_sample_carbons = len(sample)

    # Iterate through reference compound columns in the database (skipping the metadata label column)
    for col in tqdm(database.columns[1:], desc="Matching compounds"):
        db_values = pd.to_numeric(database[col], errors='coerce').dropna().values
        if len(db_values) == 0:
            continue
        
        if exact_carbon_count:
            if len(db_values) != n_sample_carbons:
                continue

        # Compute optimal shift, RMSD, and unmatched signals using dynamic programming
        rmsd, shift, unmatched = best_shift_dp(
            sample, db_values,
            shift_bounds=shift_bounds,
            grid_points=grid_points
        )

        # Store match metrics and associated metadata
        results.append({
            'compound': col,
            'RMSD': round(rmsd, 4),
            'shift': round(shift, 2),
            'solvent': database.loc[0, col],
            'carbons': len(db_values),
            'unmatched carbons': np.sort(unmatched)[::-1] if len(unmatched) else unmatched,
            'Author + Year': database.loc[1, col],
            'DOI': database.loc[2, col]
        })

    return results

import time

def run_and_compare_benchmarks(database, sample, exact_cc=False, filename="NMR_different_approaches_results"):
    """Comprehensive benchmark runner: executes all three matching methodologies (Dynamic Programming, 
    Brute Force, and Adaptive Grid Search), benchmarks their execution times, compares their top-10 ranked outputs, 
    and exports the complete results into a multi-sheet Excel workbook."""
    print("Start benchmark for the three approaches\n")
    
    # Approach 1: Monotone dynamic programming matching
    start_time = time.perf_counter()
    res_dp = compute_all_matches(sample, database, exact_carbon_count=exact_cc)
    time_dp = time.perf_counter() - start_time
    
    # Approach 2: Brute-force search
    start_time = time.perf_counter()
    res_bf = comparison_similar_compound_bruteforce(database, sample, exact_carbon_count=exact_cc)
    time_bf = time.perf_counter() - start_time
    
    # Approach 3: Adaptive grid search
    start_time = time.perf_counter()
    res_AGS = comparison_similar_compound_ags(database, sample, exact_carbon_count=exact_cc)
    time_AGS = time.perf_counter() - start_time
    
    # Convert results to DataFrames and sort by lowest RMSD error
    df_dp = pd.DataFrame(res_dp).sort_values(by="RMSD").reset_index(drop=True)
    df_bf = pd.DataFrame(res_bf).sort_values(by="RMSD").reset_index(drop=True)
    df_ags = pd.DataFrame(res_AGS).sort_values(by="RMSD").reset_index(drop=True)
    
    # Print execution time summary to console
    print("\n time: results")
    print(f"DP Monotonic Match:     {time_dp:.4f} Sekunden")
    print(f"Brute Force:            {time_bf:.4f} Sekunden")
    print(f"Adaptive Grid Search:    {time_AGS:.4f} Sekunden\n")
    
    # Compile a side-by-side comparison table of the top 10 ranked structures across all methods
    top_10_comparison = pd.DataFrame({
        "Rank": range(1, 11),
        "DP Monotonic Match": df_dp["compound"].head(10),
        "DP RMSD": df_dp["RMSD"].head(10),
        "Brute Force": df_bf["compound"].head(10),
        "BF RMSD": df_bf["RMSD"].head(10),
        "Adaptive Grid Search": df_ags["compound"].head(10),
        "AGS RMSD": df_ags["RMSD"].head(10)
    })
    
    print("Top 10 Ranking")
    print(top_10_comparison.to_string(index=False))

    # Format unmatched carbons as clean strings for Excel compatibility
    for df in [df_dp, df_bf, df_ags]:
        if 'unmatched carbons' in df.columns:
            df['unmatched carbons'] = df['unmatched carbons'].apply(lambda x: str(list(x)) if hasattr(x, '__iter__') else str(x))

    # Export all benchmark tables and top-10 comparisons into a single structured Excel workbook
    with pd.ExcelWriter(filename, engine='openpyxl') as writer:
        top_10_comparison.to_excel(writer, sheet_name="Top 10 Comparison", index=False)
        df_dp.to_excel(writer, sheet_name="DP Monotonic Match", index=False)
        df_bf.to_excel(writer, sheet_name="Brute Force", index=False)
        df_ags.to_excel(writer, sheet_name="Adaptive Grid Search", index=False)
    
    
    return df_dp, df_bf, df_ags, top_10_comparison

# how to call:
# df_dp, df_bf, df_ags, top10 = run_and_compare_benchmarks(df_database, sample_peaks, exact_cc=False, filename =)