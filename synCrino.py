import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
import numpy as np

import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import linear_sum_assignment

def load_database(file_path_database):
    database = pd.read_excel(file_path_database)
    
    df_numeric = database.iloc[3:, 1:].copy()
    
    for col in df_numeric.columns:
        cleaned_col = pd.to_numeric(df_numeric[col].astype(str).str.replace(r"\s+", "", regex=True), errors='coerce')
        sorted_values = cleaned_col.sort_values(ascending=False, na_position='last').values
        df_numeric[col] = sorted_values
        
    database.iloc[3:, 1:] = df_numeric
    
    return database

def load_sample(file_path_sample):
    sample = pd.read_excel(file_path_sample)
    sample_sorted = pd.to_numeric(sample["13 C [ppm]"], errors='coerce').round(1).dropna().sort_values(ascending=False).values
    return sample_sorted

def best_shift_differences(sample, db_values):
    sample = np.asarray(sample, dtype=float)
    db_values = np.asarray(db_values, dtype=float)

    if len(sample)==len(db_values):
        best_shift=np.average(sample-db_values)
        shifted_sample = sample - best_shift
        differences = (shifted_sample-db_values)**2
        rmsd = np.sqrt(np.mean(differences))
        unmatched_carbons = np.array([])
        
    else:
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

    
        shifts_coarse = np.linspace(-5, 5, 51) 
        for s in shifts_coarse:
            diff_matrix = ((a[:, np.newaxis] - s) - b) **2
            row_ind, col_ind = linear_sum_assignment(diff_matrix)
            score = np.sqrt(np.mean((a[row_ind] - s - b[col_ind])**2))

            if score < best_rmsd:
                best_rmsd = score
                best_shift = s

    
        shifts_medium = np.linspace(best_shift - 2*0.2, best_shift + 2*0.2, 17)
        for s in shifts_medium:
            # to make sure we're still inside the range
            if s < -5 or s > 5: continue 
        
            diff_matrix = ((a[:, np.newaxis] - s) - b) **2
            row_ind, col_ind = linear_sum_assignment(diff_matrix)
            score = np.sqrt(np.mean((a[row_ind] - s - b[col_ind])**2))

            if score < best_rmsd:
                best_rmsd = score
                best_shift = s

    
        shifts_fine = np.linspace(best_shift - 2*0.05, best_shift + 2*0.05, 21)
        for s in shifts_fine:
            if s < -5 or s > 5: continue

            diff_matrix = ((a[:, np.newaxis] - s) - b) **2
            row_ind, col_ind = linear_sum_assignment(diff_matrix)
            score = np.sqrt(np.mean((a[row_ind] - s - b[col_ind])**2))

            if score < best_rmsd:
                best_rmsd = score
                best_shift = s
            
        final_diff_matrix = ((a[:, np.newaxis] - best_shift) - b) ** 2
        a_indices, b_indices = linear_sum_assignment(final_diff_matrix)

        rmsd = np.sqrt(np.mean((a[a_indices] - best_shift - b[b_indices])**2))
  
        matched_carbons_indices = set(b_indices)
        all_b_indices = set(range(len(b)))
        unmatched_carbons_indices = list(all_b_indices - matched_carbons_indices)
        unmatched_carbons = b[unmatched_carbons_indices]
       
        

        if not sample_is_smaller:
            best_shift = -best_shift
        
    return rmsd,best_shift, unmatched_carbons

def comparison_similar_compound(database, sample, exact_carbon_count = True):
    database_sorted = database.copy()
    database_sorted = database_sorted.iloc[3:,1:]
    
    sample_sorted = np.asarray(sample, dtype=np.float64)

    if len(sample_sorted) == 0:
        raise ValueError("sample has no peaks.")
    
    n_sample_carbons = len(sample_sorted)

    results = []

    for col in database_sorted.columns:
        db_values = pd.to_numeric(database_sorted[col], errors='coerce').dropna().values
        if len(db_values) == 0:
            continue

        if exact_carbon_count:
            if len(db_values) != n_sample_carbons:
                continue
    
        rmsd, best_shift, unmatched_carbons = best_shift_differences(sample_sorted, db_values)
        
        
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
    
    results_df = pd.DataFrame(results)
    results_sorted = results_df.sort_values('RMSD', ascending=True).reset_index(drop=True)
    results_top10 = results_sorted.head(10).copy()

    results_display = results_top10.copy()
    results_display.index = results_display.index + 1

    print("Number of carbons in the sample: ", len(sample))
    print("\nTop 10 similar compounds:\n")
    print(results_display.to_string())
    

    plt.style.use('Solarize_Light2')
    fig, ax = plt.subplots(figsize=(14,4))

    state = {"idx":0}

    stem_matched = ax.stem([0], [0], label='Matched Carbons')
    stem_unmatched = ax.stem([0], [0], linefmt='r--', markerfmt='r^', basefmt=" ", label='Unmatched Carbons')
    
    stats_legend = ax.legend([Line2D([0],[0], linestyle="none")], [""], loc='upper left', bbox_to_anchor=(1.02, 0.5), frameon=True, handlelength=0)
    ax.add_artist(stats_legend)
    
    line_legend_elements = [
        Line2D([0], [0], color='#268bd2', linestyle='-', label='Matched Carbons'),
        Line2D([0], [0], color='r', linestyle='--', label='Unmatched Carbons')
    ]
    ax.legend(handles=line_legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), title="Peak types")

    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_xlim(220, 0)

    def draw(idx):
        row = results_top10.iloc[idx]
        compound = row['compound']
        unmatched_carbons = np.array(row['unmatched carbons'])
        best_shift = row['shift']

        db_series = database_sorted[compound]
        db_values = pd.to_numeric(db_series, errors='coerce').dropna().to_numpy(dtype=np.float64)
        
        shifted_sample = sample_sorted - best_shift
        final_diff_matrix = np.abs(shifted_sample[:, np.newaxis] - db_values)
        sample_indices, db_indices = linear_sum_assignment(final_diff_matrix)
        
        sample_filtered = sample_sorted[sample_indices]
        db_values_filtered = db_values[db_indices]
        
        differences = sample_filtered - db_values_filtered

        segments_matched = [np.array([[x, 0], [x, diff]]) for x, diff in zip(sample_filtered, differences)]
        stem_matched.stemlines.set_segments(segments_matched)
        
        if len(sample_filtered) > 0:
            stem_matched.markerline.set_data(sample_filtered, differences)
        else:
            stem_matched.markerline.set_data([], [])

        unmatched_height = 3
        if unmatched_carbons.size > 0:
            segments_unmatched = [np.array([[x, 0], [x, unmatched_height]]) for x in unmatched_carbons]
            stem_unmatched.stemlines.set_segments(segments_unmatched)
            stem_unmatched.markerline.set_data(unmatched_carbons, np.ones_like(unmatched_carbons) * unmatched_height)
        else:
            stem_unmatched.stemlines.set_segments([])
            stem_unmatched.markerline.set_data([], [])

        y_values = differences
        if unmatched_carbons.size > 0:
            y_values = np.concatenate([differences, np.ones_like(unmatched_carbons) * unmatched_height])
        margin = 0.5
        ax.set_ylim(np.min(y_values) - margin, np.max(y_values) + margin)

        ax.relim()
        ax.autoscale_view(scalex=False, scaley=True)

        if unmatched_carbons.size > 0:
            unmatched_formatted = [f"{val:.1f}" for val in unmatched_carbons]
            unmatched_str = f"[{', '.join(unmatched_formatted)}]"
        else:
            unmatched_str = "[]"
                
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
        if event.key == "right":
            state["idx"] = min(state["idx"] + 1, len(results_top10) - 1)
            draw(state["idx"])
        elif event.key == "left":
            state["idx"] = max(state["idx"] - 1, 0)
            draw(state["idx"])
        #elif event.key == "s":  # Drücke 's' auf der Tastatur zum Speichern des aktuellen Rangs
            #idx = state["idx"]
            #compound_name = results_top10.iloc[idx]['compound']
            #filename = f"Match_Rank_{idx+1}_{compound_name}.tiff"
            #plt.savefig(filename, dpi=600)
            #print(f"Saved: {filename}")

    fig.canvas.mpl_connect("key_press_event", on_key)

    draw(0)

    plt.tight_layout()
    plt.subplots_adjust(right=0.65)
    #plt.savefig('sample_hyperforin_vs_hyperforin.tiff', dpi =600)
    
    plt.show()

    return results_df
