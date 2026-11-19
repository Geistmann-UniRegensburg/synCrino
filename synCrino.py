# synCrino: 13C NMR chemical shift matching tool
#Description: Automated alignment and comparison of experimental 13C NMR chemical shifts against reference database entries.

import pandas as pd
# Suppress pandas warning
pd.set_option('future.no_silent_downcasting', True)
import numpy as np

import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import linear_sum_assignment


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

def best_shift_differences(sample, db_values):
    """Find the optimal chemical shift correction, RMSD, and unmatched signals between sample and database sets."""
    sample = np.asarray(sample, dtype=float)
    db_values = np.asarray(db_values, dtype=float)

    if len(sample)==len(db_values):
        # For equal-length sets (same carbon count), calculate the uniform systematic shift correction via the average difference
        best_shift=np.average(sample-db_values)
        shifted_sample = sample - best_shift

        # Compute RMSD (Root Mean Square Deviation) as the primary metric for database ranking
        differences = (shifted_sample-db_values)**2
        rmsd = np.sqrt(np.mean(differences))

        # No unmatched signals exist when carbon counts are identical
        unmatched_carbons = np.array([])
        
    else:
        # For non-equal-length sets (different carbon counts)
        # Assign 'a' to the smaller array and 'b' to the larger array for optimization
        if(len(sample)) < len(db_values):
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
        # optimal bipartite matching problem. We use the Hungarian algorithm to pair 
        # signals and minimize the RMSD across a hierarchical 3-step grid search (coarse -> medium -> fine resolution):

        # 1. Coarse grid search: Test a wide range of potential systematic shifts (-5 to +5 ppm) 
        shifts_coarse = np.linspace(-5, 5, 51)
        for s in shifts_coarse:
            diff_matrix = ((a[:, np.newaxis] - s) - b) **2
            row_ind, col_ind = linear_sum_assignment(diff_matrix)
            score = np.sqrt(np.mean((a[row_ind] - s - b[col_ind])**2))

            if score < best_rmsd:
                best_rmsd = score
                best_shift = s

        # 2. Medium resolution grid search: Refine search window around the best coarse shift
        shifts_medium = np.linspace(best_shift - 2*0.2, best_shift + 2*0.2, 17)
        for s in shifts_medium:
            # Ensure the evaluated shift stays within valid boundaries
            if s < -5 or s > 5: continue 
        
            diff_matrix = ((a[:, np.newaxis] - s) - b) **2
            row_ind, col_ind = linear_sum_assignment(diff_matrix)
            score = np.sqrt(np.mean((a[row_ind] - s - b[col_ind])**2))

            if score < best_rmsd:
                best_rmsd = score
                best_shift = s

        # 3. Fine resolution grid search: Final high-precision optimization around the best medium shift
        shifts_fine = np.linspace(best_shift - 2*0.05, best_shift + 2*0.05, 21)
        for s in shifts_fine:
            if s < -5 or s > 5: continue

            diff_matrix = ((a[:, np.newaxis] - s) - b) **2
            row_ind, col_ind = linear_sum_assignment(diff_matrix)
            score = np.sqrt(np.mean((a[row_ind] - s - b[col_ind])**2))

            if score < best_rmsd:
                best_rmsd = score
                best_shift = s

        # Final assignment using the optimized shift    
        final_diff_matrix = ((a[:, np.newaxis] - best_shift) - b) ** 2
        a_indices, b_indices = linear_sum_assignment(final_diff_matrix)

        # Final calculation of RMSD
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

def comparison_similar_compound(database, sample, exact_carbon_count = True):
    """Compare experimental sample 13C NMR shifts against a reference database, rank by RMSD, and provide interactive visualization.
    By default it only compares the sample with database entries with the same amount of carbons; it can be switch by setting it to False"""
    database_sorted = database.copy()
    database_sorted = database_sorted.iloc[3:,1:]
    
    sample_sorted = np.asarray(sample, dtype=np.float64)

    # Ensure sample has values
    if len(sample_sorted) == 0:
        raise ValueError("sample has no peaks.")
    
    n_sample_carbons = len(sample_sorted)

    results = []

    # Iterate through database columns (compounds) to compute match scores
    for col in database_sorted.columns:
        db_values = pd.to_numeric(database_sorted[col], errors='coerce').dropna().values
        if len(db_values) == 0:
            continue

        # Optionally filter out compounds that don't match the exact carbon count
        if exact_carbon_count:
            if len(db_values) != n_sample_carbons:
                continue

        # Compute optimal shift correction, RMSD, and unmatched signals
        rmsd, best_shift, unmatched_carbons = best_shift_differences(sample_sorted, db_values)
        
        # Store comprehensive match metrics and metadata for ranking
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

    # Compile results and extract the top 10 best matches ranked by RMSD
    results_df = pd.DataFrame(results)
    results_sorted = results_df.sort_values('RMSD', ascending=True).reset_index(drop=True)
    results_top10 = results_sorted.head(10).copy()

    results_display = results_top10.copy()
    results_display.index = results_display.index + 1

    print("Number of carbons in the sample: ", len(sample))
    print("\nTop 10 similar compounds:\n")
    print(results_display.to_string())
    
    # Initialize interactive Matplotlib visualization for the top matches
    plt.style.use('Solarize_Light2')
    fig, ax = plt.subplots(figsize=(14,4))

    state = {"idx":0}

    # Initialize stem plots for matched and unmatched carbon signals
    stem_matched = ax.stem([0], [0], label='Matched Carbons')
    stem_unmatched = ax.stem([0], [0], linefmt='r--', markerfmt='r^', basefmt=" ", label='Unmatched Carbons')
    
    # Set up dedicated text box legend for compound metadata statistics
    stats_legend = ax.legend([Line2D([0],[0], linestyle="none")], [""], loc='upper left', bbox_to_anchor=(1.02, 0.5), frameon=True, handlelength=0)
    ax.add_artist(stats_legend)
    
    line_legend_elements = [
        Line2D([0], [0], color='#268bd2', linestyle='-', label='Matched Carbons'),
        Line2D([0], [0], color='r', linestyle='--', label='Unmatched Carbons')
    ]
    ax.legend(handles=line_legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), title="Peak types")

    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_xlim(220, 0) # Standard NMR chemical shift axis convention (high to low ppm)

    def draw(idx):
        """Update plot elements for the currently selected match rank."""
        row = results_top10.iloc[idx]
        compound = row['compound']
        unmatched_carbons = np.array(row['unmatched carbons'])
        best_shift = row['shift']

        # Retrieve reference carbon shifts for the selected compound
        db_series = database_sorted[compound]
        db_values = pd.to_numeric(db_series, errors='coerce').dropna().to_numpy(dtype=np.float64)
        
        # Apply optimal shift and re-run assignment for plotting residuals
        shifted_sample = sample_sorted - best_shift
        final_diff_matrix = np.abs(shifted_sample[:, np.newaxis] - db_values)
        sample_indices, db_indices = linear_sum_assignment(final_diff_matrix)
        
        sample_filtered = sample_sorted[sample_indices]
        db_values_filtered = db_values[db_indices]
        
        differences = sample_filtered - db_values_filtered

        # Update matched stems visualization segments
        segments_matched = [np.array([[x, 0], [x, diff]]) for x, diff in zip(sample_filtered, differences)]
        stem_matched.stemlines.set_segments(segments_matched)
        
        if len(sample_filtered) > 0:
            stem_matched.markerline.set_data(sample_filtered, differences)
        else:
            stem_matched.markerline.set_data([], [])

        # Update unmatched stems visualization segments (rendered at a fixed height)
        unmatched_height = 3
        if unmatched_carbons.size > 0:
            segments_unmatched = [np.array([[x, 0], [x, unmatched_height]]) for x in unmatched_carbons]
            stem_unmatched.stemlines.set_segments(segments_unmatched)
            stem_unmatched.markerline.set_data(unmatched_carbons, np.ones_like(unmatched_carbons) * unmatched_height)
        else:
            stem_unmatched.stemlines.set_segments([])
            stem_unmatched.markerline.set_data([], [])

        # Dynamically adjust Y-axis limits based on residuals and unmatched markers
        y_values = differences
        if unmatched_carbons.size > 0:
            y_values = np.concatenate([differences, np.ones_like(unmatched_carbons) * unmatched_height])
        margin = 0.5
        ax.set_ylim(np.min(y_values) - margin, np.max(y_values) + margin)

        ax.relim()
        ax.autoscale_view(scalex=False, scaley=True)

        # Format unmatched carbons string for metadata display
        if unmatched_carbons.size > 0:
            unmatched_formatted = [f"{val:.1f}" for val in unmatched_carbons]
            unmatched_str = f"[{', '.join(unmatched_formatted)}]"
        else:
            unmatched_str = "[]"
                
        # Update metadata text box content
        stats_text = (
            f"compound: {compound}\n"
            f"RMSD: {row['RMSD']:.4f}\n"
            f"shift: {row['shift']:.2f}\n"
            f"solvent: {row['solvent']}\n"
            f"Cs sample vs Cs compound: {len(sample_sorted)} vs {row['carbons']}\n"
            f"unmatched carbons: {unmatched_str}"
        )
        stats_legend.get_texts()[0].set_text(stats_text)

        ax.set_title(f"Match Rank {idx+1} (of Top 10)")
        fig.canvas.draw_idle()

    def on_key(event):
        """Handle keyboard arrow events to cycle through top matches."""
        if event.key == "right":
            state["idx"] = min(state["idx"] + 1, len(results_top10) - 1)
            draw(state["idx"])
        elif event.key == "left":
            state["idx"] = max(state["idx"] - 1, 0)
            draw(state["idx"])
        

    fig.canvas.mpl_connect("key_press_event", on_key)

    # Render the first rank by default
    draw(0)

    plt.tight_layout()
    plt.subplots_adjust(right=0.65)
    
    plt.show()

    return results_df
