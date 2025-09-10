import streamlit as st
import pandas as pd
import numpy as np
import os
import glob
import tempfile
import zipfile
from pathlib import Path
import re
import shutil
from io import BytesIO
import difflib

# Set page configuration
st.set_page_config(
    page_title="Ct Value Extraction Tool",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --------------------------
# DATA PROCESSING FUNCTIONS
# --------------------------

def process_single_csv(file_path):
    """Process a single CSV file and extract sample and positive control data"""
    try:
        # Read the first two rows to get headers
        header_rows = pd.read_csv(file_path, header=None, nrows=2)
        
        # Extract headers
        header1 = header_rows.iloc[0].fillna('').astype(str).tolist()
        header2 = header_rows.iloc[1].fillna('').astype(str).tolist()
        
        # Combine headers
        full_header = header1.copy()
        for i in range(5, min(21, len(header1))):
            if i < len(header2) and header2[i]:
                full_header[i] = header2[i]
        
        # Read the main data
        raw_data = pd.read_csv(file_path, header=None, skiprows=2)
        
        # Set column names
        max_cols = min(len(full_header), len(raw_data.columns))
        raw_data.columns = [f"col_{i}" for i in range(len(raw_data.columns))]  # Temporary column names
        for i in range(max_cols):
            if i < len(full_header):
                raw_data.rename(columns={f"col_{i}": full_header[i]}, inplace=True)
        
        # Identify pathogen columns
        pathogen_cols = list(range(5, min(21, len(raw_data.columns)), 2)
        pathogen_names = [full_header[i] for i in pathogen_cols if i < len(full_header)]
        
        # Create result and Ct column names
        result_ct_names = []
        for pathogen in pathogen_names:
            result_ct_names.extend([f"{pathogen}_Result", f"{pathogen}_Ct"])
        
        # Set all column names
        all_cols = full_header[:5] + result_ct_names
        if len(raw_data.columns) > 21:
            all_cols.extend(full_header[21:len(raw_data.columns)])
        
        # Ensure we have enough column names
        if len(all_cols) < len(raw_data.columns):
            all_cols.extend([f"col_{i}" for i in range(len(all_cols), len(raw_data.columns))])
        
        raw_data.columns = all_cols[:len(raw_data.columns)]
        
        # Add source file information
        raw_data['Source_File'] = os.path.basename(file_path)
        
        # Process positive controls
        pc_mask = (raw_data.get('Type', '') == 'PC') | (
            'Auto   Interpretation' in raw_data.columns and 
            raw_data['Auto   Interpretation'].ast(str).str.contains('Positive Control', na=False)
        )
        
        pc_data = raw_data[pc_mask].copy()
        
        # Process PC data
        pc_results = []
        for pathogen in pathogen_names:
            result_col = f"{pathogen}_Result"
            ct_col = f"{pathogen}_Ct"
            
            if result_col in pc_data.columns and ct_col in pc_data.columns:
                positive_pcs = pc_data[pc_data[result_col] == '+']
                if not positive_pcs.empty:
                    pc_ct = positive_pcs[ct_col].iloc[0]
                    pc_status = "Passed" if not pd.isna(pc_ct) else "Failed"
                    pc_results.append({
                        'Pathogen': pathogen,
                        'PC_Ct': pc_ct,
                        'PC_Status': pc_status,
                        'Source_File': os.path.basename(file_path)
                    })
        
        pc_df = pd.DataFrame(pc_results) if pc_results else pd.DataFrame()
        
        return {
            'sample_data': raw_data,
            'pc_data': pc_df,
            'pc_status': all([r['PC_Status'] == 'Passed' for r in pc_results]) if pc_results else False,
            'all_pathogens': pathogen_names
        }
        
    except Exception as e:
        st.error(f"Error processing file {file_path}: {str(e)}")
        return None

def extract_zip_and_process(uploaded_zip, selected_samples):
    """Extract zip file and process all CSV files"""
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Save the uploaded zip file
        zip_path = os.path.join(tmp_dir, "uploaded_files.zip")
        with open(zip_path, "wb") as f:
            f.write(uploaded_zip.getvalue())
        
        # Extract the zip file
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(tmp_dir)
        
        # Find all CSV files
        csv_files = []
        for root, dirs, files in os.walk(tmp_dir):
            for file in files:
                if file.endswith('.csv'):
                    csv_files.append(os.path.join(root, file))
        
        if not csv_files:
            raise ValueError("No CSV files found in the uploaded zip.")
        
        all_sample_data = []
        all_pc_data = []
        pc_status_list = []
        not_found_files = []
        
        for file_path in csv_files:
            processed = process_single_csv(file_path)
            if processed is not None:
                all_sample_data.append(processed['sample_data'])
                if not processed['pc_data'].empty:
                    all_pc_data.append(processed['pc_data'])
                pc_status_list.append({
                    'File': os.path.basename(file_path),
                    'PC_Worked': processed['pc_status'],
                    'Message': "All controls passed" if processed['pc_status'] else "Some controls failed"
                })
            else:
                not_found_files.append(os.path.basename(file_path))
        
        # Combine all data
        combined_samples = pd.concat(all_sample_data, ignore_index=True) if all_sample_data else pd.DataFrame()
        combined_pc = pd.concat(all_pc_data, ignore_index=True) if all_pc_data else pd.DataFrame()
        pc_summary = pd.DataFrame(pc_status_list) if pc_status_list else pd.DataFrame()
        
        # Clean sample names
        if 'Name' in combined_samples.columns:
            combined_samples['Name'] = combined_samples['Name'].astype(str).str.strip()
        
        # Match samples
        filtered_samples = match_samples(combined_samples, selected_samples)
        
        # Select relevant columns
        ct_cols = [col for col in filtered_samples.columns if col.endswith('_Ct')] if not filtered_samples.empty else []
        selected_cols = ['Name'] + ct_cols
        if 'Auto   Interpretation' in filtered_samples.columns:
            selected_cols.append('Auto   Interpretation')
        
        # Ensure Source_File is included
        if 'Source_File' not in selected_cols and 'Source_File' in filtered_samples.columns:
            selected_cols = ['Source_File'] + selected_cols
        
        # Filter to only include available columns
        available_cols = [col for col in selected_cols if col in filtered_samples.columns]
        filtered_samples = filtered_samples[available_cols] if available_cols else pd.DataFrame()
        
        return {
            'sample_results': filtered_samples,
            'pc_details': combined_pc,
            'pc_summary': pc_summary,
            'all_pathogens': combined_pc['Pathogen'].unique().tolist() if not combined_pc.empty else [],
            'not_found_files': not_found_files,
            'csv_count': len(csv_files)
        }

def match_samples(combined_samples, selected_samples):
    """Match samples using Python's built-in difflib for fuzzy matching"""
    if combined_samples.empty or 'Name' not in combined_samples.columns:
        return pd.DataFrame()
    
    # Clean names
    combined_samples['Name_clean'] = combined_samples['Name'].astype(str).str.upper().str.strip()
    selected_clean = [str(s).upper().strip() for s in selected_samples]
    
    # Get all unique clean names from the data
    all_names = combined_samples['Name_clean'].unique().tolist()
    
    # Perform fuzzy matching using difflib
    matched_samples = []
    for episode in selected_clean:
        # Find the best match using difflib
        matches = difflib.get_close_matches(episode, all_names, n=1, cutoff=0.8)
        
        if matches:
            best_match = matches[0]
            matched_data = combined_samples[combined_samples['Name_clean'] == best_match].copy()
            matched_data['Episode_Number'] = episode
            matched_samples.append(matched_data)
    
    if matched_samples:
        return pd.concat(matched_samples, ignore_index=True)
    else:
        return pd.DataFrame()

# --------------------------
# STREAMLIT UI
# --------------------------

# Initialize session state
if 'results' not in st.session_state:
    st.session_state.results = None
if 'pc_summary' not in st.session_state:
    st.session_state.pc_summary = None
if 'pc_details' not in st.session_state:
    st.session_state.pc_details = None
if 'not_found_files' not in st.session_state:
    st.session_state.not_found_files = []
if 'processed' not in st.session_state:
    st.session_state.processed = False
if 'csv_count' not in st.session_state:
    st.session_state.csv_count = 0

# Sidebar
with st.sidebar:
    st.title("Ct Value Extraction Tool")
    
    st.header("Data Input")
    
    # Use file uploader for CSV files (zip format)
    st.info("Upload a ZIP file containing all CSV files")
    uploaded_zip = st.file_uploader(
        "Upload ZIP file with CSV files:",
        type=["zip"],
        help="Create a ZIP file containing all your CSV files and upload it here"
    )
    
    # Sample list upload
    sample_file = st.file_uploader(
        "Upload Sample List (with 'Episode Number' column):",
        type=["csv", "xlsx", "xls", "txt"]
    )
    
    st.header("Processing Options")
    include_pc = st.checkbox("Include Positive Control Results in Download", value=True)
    file_format = st.radio("Download Format:", options=["CSV", "Excel"], index=0)
    
    process_btn = st.button("Extract Ct Values", type="primary")
    
    st.header("Download")

# Main content
st.title("Ct Value Extraction Tool")

# Process files when button is clicked
if process_btn and uploaded_zip and sample_file:
    with st.spinner("Processing files..."):
        # Read sample list
        try:
            if sample_file.name.endswith('.csv'):
                sample_df = pd.read_csv(sample_file)
            elif sample_file.name.endswith(('.xlsx', '.xls')):
                sample_df = pd.read_excel(sample_file)
            else:
                sample_df = pd.read_csv(sample_file, delimiter='\t')
            
            if 'Episode Number' not in sample_df.columns:
                st.error("Sample file must contain an 'Episode Number' column")
                st.stop()
            
            selected_samples = sample_df['Episode Number'].tolist()
            
            # Process files
            results = extract_zip_and_process(uploaded_zip, selected_samples)
            
            # Store results in session state
            st.session_state.results = results['sample_results']
            st.session_state.pc_summary = results['pc_summary']
            st.session_state.pc_details = results['pc_details']
            st.session_state.not_found_files = results['not_found_files']
            st.session_state.csv_count = results['csv_count']
            st.session_state.processed = True
            
            st.success("Processing complete!")
            
        except Exception as e:
            st.error(f"Error processing files: {str(e)}")

# Display status information
st.header("Processing Status")

if st.session_state.processed:
    n_processed = len(st.session_state.results) if st.session_state.results is not None else 0
    
    if sample_file:
        try:
            if sample_file.name.endswith('.csv'):
                sample_df = pd.read_csv(sample_file)
            elif sample_file.name.endswith(('.xlsx', '.xls')):
                sample_df = pd.read_excel(sample_file)
            else:
                sample_df = pd.read_csv(sample_file, delimiter='\t')
            
            n_selected = len(sample_df)
            
            # Check for missing samples
            if st.session_state.results is not None and 'Episode_Number' in st.session_state.results.columns:
                processed_episodes = st.session_state.results['Episode_Number'].unique().tolist()
                selected_episodes = sample_df['Episode Number'].astype(str).str.upper().str.strip().tolist()
                missing_samples = set(selected_episodes) - set(processed_episodes)
            else:
                missing_samples = set()
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Selected Samples", n_selected)
            with col2:
                st.metric("CSV Files Processed", st.session_state.csv_count)
            with col3:
                st.metric("Matched Samples", n_processed)
            
            if missing_samples:
                st.warning(f"⚠️ {len(missing_samples)} samples not found: {', '.join(list(missing_samples)[:5])}{'...' if len(missing_samples) > 5 else ''}")
            
            if st.session_state.not_found_files:
                st.warning(f"⚠️ {len(st.session_state.not_found_files)} files could not be processed")
        
        except:
            pass
else:
    st.info("Upload a ZIP file with CSV files and a sample list, then click 'Extract Ct Values' to process.")

# Create tabs
tab1, tab2, tab3 = st.tabs(["Sample Preview", "Positive Controls", "Download"])

# Sample Preview tab
with tab1:
    st.header("Sample Results Preview")
    
    if st.session_state.results is not None and not st.session_state.results.empty:
        # Display the data
        st.dataframe(
            st.session_state.results,
            use_container_width=True,
            hide_index=True
        )
    else:
        st.info("No sample data to display. Process files first.")

# Positive Controls tab
with tab2:
    st.header("Positive Control Summary")
    
    if st.session_state.pc_summary is not None and not st.session_state.pc_summary.empty:
        st.dataframe(
            st.session_state.pc_summary,
            use_container_width=True,
            hide_index=True
        )
    else:
        st.info("No positive control data to display.")
    
    st.header("Detailed Positive Control Results")
    
    if st.session_state.pc_details is not None and not st.session_state.pc_details.empty:
        st.dataframe(
            st.session_state.pc_details,
            use_container_width=True,
            hide_index=True
        )
    else:
        st.info("No detailed positive control data to display.")

# Download tab
with tab3:
    st.header("Download Results")
    
    if st.session_state.processed and st.session_state.results is not None:
        # Prepare download data
        if file_format == "Excel":
            # Create Excel file
            output = BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                st.session_state.results.to_excel(writer, sheet_name='Sample_Results', index=False)
                if include_pc:
                    if st.session_state.pc_summary is not None and not st.session_state.pc_summary.empty:
                        st.session_state.pc_summary.to_excel(writer, sheet_name='PC_Summary', index=False)
                    if st.session_state.pc_details is not None and not st.session_state.pc_details.empty:
                        st.session_state.pc_details.to_excel(writer, sheet_name='PC_Details', index=False)
            
            # Create download button
            st.download_button(
                label="Download Excel File",
                data=output.getvalue(),
                file_name=f"Ct_results_{pd.Timestamp.now().strftime('%Y-%m-%d')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        else:
            # Create CSV files and zip them
            zip_buffer = BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                # Add sample results
                sample_csv = st.session_state.results.to_csv(index=False)
                zip_file.writestr(f"Ct_results_{pd.Timestamp.now().strftime('%Y-%m-%d')}.csv", sample_csv)
                
                # Add PC data if requested
                if include_pc:
                    if st.session_state.pc_summary is not None and not st.session_state.pc_summary.empty:
                        pc_summary_csv = st.session_state.pc_summary.to_csv(index=False)
                        zip_file.writestr(f"PC_summary_{pd.Timestamp.now().strftime('%Y-%m-%d')}.csv", pc_summary_csv)
                    
                    if st.session_state.pc_details is not None and not st.session_state.pc_details.empty:
                        pc_details_csv = st.session_state.pc_details.to_csv(index=False)
                        zip_file.writestr(f"PC_details_{pd.Timestamp.now().strftime('%Y-%m-%d')}.csv", pc_details_csv)
            
            # Create download button
            st.download_button(
                label="Download ZIP File",
                data=zip_buffer.getvalue(),
                file_name=f"Ct_results_{pd.Timestamp.now().strftime('%Y-%m-%d')}.zip",
                mime="application/zip"
            )
    else:
        st.info("No data available for download. Process files first.")
