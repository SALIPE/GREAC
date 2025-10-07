import json
import os
import shutil
import subprocess
import zipfile
from pathlib import Path
from typing import Dict, List, Any
from glob import glob

import streamlit as st

st.set_page_config(
    page_title="GREAC",
    page_icon="🧬",
    layout="wide"
)

# ============================================================================
# SCRIPTS CONFIGURATION
# ============================================================================
GROUPS_BASE_PATH = "~/.project_cache"
SCRIPTS_CONFIG = {
    "Benchmark GREAC": {
        "script_path": "../scripts/local/benchmark.sh",
        "description": "Executes the complete GREAC benchmark",
        "parameters": [
            {
                "name": "train_dir",
                "label": "📂 Training Directory",
                "type": "text",
                "default": "~/Desktop/datasets/dengue/train/kmers",
                "required": True
            },
            {
                "name": "test_dir",
                "label": "📂 Test Directory",
                "type": "text",
                "default": "~/Desktop/datasets/dengue/test",
                "required": True
            },
            {
                "name": "group_name",
                "label": "🏷 Group Name",
                "type": "group_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "window_size",
                "label": "🪟 Window Size",
                "type": "window_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "metric",
                "label": "📊 Metric",
                "type": "select",
                "options": ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"],
                "default": "manhattan",
                "required": True
            },
            {
                "name": "kmer_size",
                "label": "🔢 K-mer Size",
                "type": "number",
                "default": 7,
                "step": 1,
                "required": True
            },
            {
                "name": "threshold",
                "label": "📏 Window Consideration Percentage",
                "type": "number",
                "default": 0.0,
                "step": 0.1,
                "format": "%.1f",
                "required": True
            },
            {
                "name": "reference_path",
                "label": "📄 Reference File (.fasta)",
                "type": "file_path",
                "extensions": [".fasta", ".fa", ".fas"],
                "default": "~/Desktop/datasets/denv/refseq.fasta",
                "required": True
            }
        ]
    },
    "FASTA Regions": {
        "script_path": "../scripts/local/fasta-regions.sh",
        "description": "Returns the .bed file with regions and reduces the dataset (if provided)",
        "parameters": [
            {
                "name": "group_name",
                "label": "🏷 Group Name",
                "type": "group_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "window_size",
                "label": "🪟 Window Size",
                "type": "window_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "dataset",
                "label": "📂 Input Dataset Directory (Optional)",
                "type": "text",
                "default": "",
                "required": False,
                "help": "Leave empty if you don't want to reduce a dataset"
            }
        ]
    },
    "Feature Extraction": {
        "script_path": "../scripts/local/extract-features.sh",
        "description": "Applies the GREAC method for region extraction and model training",
        "parameters": [
            {
                "name": "train_dir",
                "label": "📂 Training Directory",
                "type": "text",
                "default": "~/Desktop/datasets/dengue/train/kmers",
                "required": True
            },
            {
                "name": "group_name",
                "label": "🏷 Group Name",
                "type": "group_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "window_size",
                "label": "🪟 Window Size",
                "type": "window_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "metric",
                "label": "📊 Metric",
                "type": "select",
                "options": ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"],
                "default": "manhattan",
                "required": True
            },
            {
                "name": "threshold",
                "label": "📏 Window Consideration Percentage",
                "type": "number",
                "default": 0.0,
                "step": 0.1,
                "format": "%.1f",
                "required": True
            },
            {
                "name": "reference_path",
                "label": "📄 Reference File (.fasta)",
                "type": "file_path",
                "extensions": [".fasta", ".fa", ".fas"],
                "default": "~/Desktop/datasets/denv/refseq.fasta",
                "required": True
            }
        ]
    },
    "Classification": {
        "script_path": "../scripts/local/file-classification.sh",
        "description": "Classifies a file based on a pre-trained model",
        "parameters": [
            {
                "name": "train_dir",
                "label": "📂 Training Directory",
                "type": "text",
                "default": "~/Desktop/datasets/dengue/train/kmers",
                "required": True
            },
            {
                "name": "test_file",
                "label": "📄 Test File",
                "type": "file_path",
                "extensions": [".fasta", ".fa", ".fas"],
                "default": "~/Desktop/datasets/dengue/test/sample.fasta",
                "required": True,
                "help": "Single FASTA file to classify"
            },
            {
                "name": "group_name",
                "label": "🏷 Group Name",
                "type": "group_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "window_size",
                "label": "🪟 Window Size",
                "type": "window_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "metric",
                "label": "📊 Metric",
                "type": "select",
                "options": ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"],
                "default": "manhattan",
                "required": True
            },
            {
                "name": "kmer_size",
                "label": "🔢 K-mer Size",
                "type": "number",
                "default": 7,
                "step": 1,
                "required": True
            },
            {
                "name": "threshold",
                "label": "📏 Window Consideration Percentage",
                "type": "number",
                "default": 0.0,
                "step": 0.1,
                "format": "%.1f",
                "required": True
            },
            {
                "name": "reference_path",
                "label": "📄 Reference File (.fasta)",
                "type": "file_path",
                "extensions": [".fasta", ".fa", ".fas"],
                "default": "~/Desktop/datasets/denv/refseq.fasta",
                "required": True
            }
        ]
    },
    "Complete Pipeline (doall.sh)": {
        "script_path": "../scripts/local/doall.sh",
        "description": "Executes the complete GREAC pipeline: data preparation, k-mer extraction, and benchmarking",
        "parameters": [
            {
                "name": "group_name",
                "label": "🏷 Group Name",
                "type": "select",
                "options": ["denv", "hbv", "hiv", "hiv2", "sars", "monkeypox", 
                           "bees1", "bees2", "bees3", "bees4", "bees5", "bees6", 
                           "bees7", "bees8", "bees9", "bees10", "bees11", "bees12",
                           "bees13", "bees14", "bees15", "bees16"],
                "default": "denv",
                "required": True,
                "help": "Dataset to process (predefined paths in script)"
            },
            {
                "name": "window_size",
                "label": "🪟 Window Size",
                "type": "text",
                "default": "500",
                "required": True
            },
            {
                "name": "kmer_size",
                "label": "🔢 K-mer Size",
                "type": "number",
                "default": 7,
                "step": 1,
                "required": True
            },
            {
                "name": "threshold",
                "label": "📏 Window Consideration Percentage",
                "type": "number",
                "default": 0.0,
                "step": 0.1,
                "format": "%.1f",
                "required": True
            }
        ]
    }
}

# Directories/patterns where outputs can be found (supports wildcards)
OUTPUT_DIRECTORIES = [
    "~/.project_cache",
    "../GREAC/output-*",
    "../GREAC/*.pdf",
]

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def get_available_groups(base_path: str) -> List[str]:
    """Lists available directories as groups"""
    expanded_path = os.path.expanduser(base_path)
    if not os.path.exists(expanded_path):
        return []
    
    try:
        groups = [d for d in os.listdir(expanded_path) 
                 if os.path.isdir(os.path.join(expanded_path, d)) and not d.startswith('.')]
        return sorted(groups)
    except:
        return []

def get_available_windows(base_path: str, group_name: str) -> List[str]:
    """Lists available subdirectories as windows within a group"""
    if not group_name:
        return []
    
    expanded_path = os.path.expanduser(base_path)
    if not os.path.exists(expanded_path):
        return []
    
    group_path = os.path.join(expanded_path, group_name)
    if not os.path.exists(group_path):
        return []
    
    try:
        windows = [d for d in os.listdir(group_path) 
                  if os.path.isdir(os.path.join(group_path, d)) and not d.startswith('.')]
        return sorted(windows)
    except:
        return []
    
def list_directory_tree(directory: str, prefix: str = "", max_depth: int = 5) -> str:
    """Lists the directory structure in tree format"""
    expanded_dir = os.path.expanduser(directory)
    if not os.path.exists(expanded_dir):
        return f"❌ Directory not found: {directory}\n"
    
    tree_str = f"📁 {os.path.basename(expanded_dir)}/\n"
    
    try:
        items = []
        for root, dirs, files in os.walk(expanded_dir):
            level = root.replace(expanded_dir, "").count(os.sep)
            if level >= max_depth:
                continue
                
            indent = "    " * level
            
            # Add directories
            for d in sorted(dirs):
                items.append(f"{indent}├── 📁 {d}/")
            
            # Add files
            for f in sorted(files):
                file_path = os.path.join(root, f)
                size = os.path.getsize(file_path)
                size_str = format_file_size(size)
                items.append(f"{indent}├── 📄 {f} ({size_str})")
        
        tree_str += "\n".join(items)
    except Exception as e:
        tree_str += f"❌ Error listing: {str(e)}"
    
    return tree_str

def format_file_size(size_bytes: int) -> str:
    """Formats file size in readable format"""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"

def render_parameter(param_config: Dict[str, Any], key_prefix: str, context: Dict[str, Any] = None) -> Any:
    """Renders a parameter according to its type"""
    param_type = param_config.get("type", "text")
    name = param_config["name"]
    label = param_config.get("label", name)
    default = param_config.get("default", "")
    required = param_config.get("required", False)
    help_text = param_config.get("help", None)
    
    key = f"{key_prefix}_{name}"
    
    if param_type == "text":
        value = st.text_input(
            label + (" *" if required else ""),
            value=default,
            key=key,
            help=help_text
        )
    elif param_type == "group_selector":
        # Hybrid field for group selection/typing
        base_path = os.path.expanduser(GROUPS_BASE_PATH)
        available_groups = get_available_groups(base_path)
        
        st.markdown(f"**{label}**" + (" *" if required else ""))
        if help_text:
            st.caption(help_text)
        
        col1, col2 = st.columns([3, 1])
        with col1:
            if available_groups:
                # Radio button to choose mode
                mode = st.radio(
                    "Input mode:",
                    options=["Select existing", "Type new"],
                    key=f"{key}_mode",
                    horizontal=True
                )
                
                if mode == "Select existing":
                    value = st.selectbox(
                        "Choose a group:",
                        options=[""] + available_groups,
                        key=f"{key}_select"
                    )
                else:
                    value = st.text_input(
                        "Enter group name:",
                        value=default,
                        key=f"{key}_text"
                    )
            else:
                st.info(f"ℹ️ No groups found in: {base_path}")
                value = st.text_input(
                    "Enter group name:",
                    value=default,
                    key=f"{key}_text"
                )
        
        with col2:
            if available_groups:
                st.metric("Available groups", len(available_groups))
    
    elif param_type == "window_selector":
        # Hybrid field for window selection/typing (depends on group)
        base_path = os.path.expanduser(GROUPS_BASE_PATH)
        group_name = context.get("group_name") if context else None
        
        st.markdown(f"**{label}**" + (" *" if required else ""))
        if help_text:
            st.caption(help_text)
        
        if not group_name:
            st.warning("⚠️ Select a group first to see available windows")
            value = st.text_input(
                "Enter value:",
                value=default,
                key=f"{key}_text"
            )
        else:
            available_windows = get_available_windows(base_path, group_name)
            
            col1, col2 = st.columns([3, 1])
            with col1:
                if available_windows:
                    # Radio button to choose mode
                    mode = st.radio(
                        "Input mode:",
                        options=["Select existing", "Type new"],
                        key=f"{key}_mode",
                        horizontal=True
                    )
                    
                    if mode == "Select existing":
                        value = st.selectbox(
                            f"Choose a window (from group '{group_name}'):",
                            options=[""] + available_windows,
                            key=f"{key}_select"
                        )
                    else:
                        value = st.text_input(
                            "Enter value:",
                            value=default,
                            key=f"{key}_text"
                        )
                else:
                    st.info(f"ℹ️ No windows found for group '{group_name}'")
                    value = st.text_input(
                        "Enter value:",
                        value=default,
                        key=f"{key}_text"
                    )
            
            with col2:
                if available_windows:
                    st.metric("Available windows", len(available_windows))
    
    elif param_type == "number":
        step = param_config.get("step", 1)
        format_str = param_config.get("format", "%d")
        value = st.number_input(
            label + (" *" if required else ""),
            value=default,
            step=step,
            format=format_str,
            key=key,
            help=help_text
        )
    elif param_type == "select":
        options = param_config.get("options", [])
        default_index = 0
        if default in options:
            default_index = options.index(default)
        value = st.selectbox(
            label + (" *" if required else ""),
            options=options,
            index=default_index,
            key=key,
            help=help_text
        )
    elif param_type == "checkbox":
        value = st.checkbox(
            label,
            value=default,
            key=key,
            help=help_text
        )
    elif param_type == "file_path":
        value = st.text_input(
            label + (" *" if required else ""),
            value=default,
            key=key,
            help=help_text or f"Accepted extensions: {', '.join(param_config.get('extensions', []))}"
        )
    else:
        value = st.text_input(label, value=default, key=key, help=help_text)
    
    return value

def build_command(script_name: str, script_path: str, parameters: List[Dict], param_values: Dict) -> List[str]:
    """Builds the command to be executed based on script type"""
    cmd = [script_path]
    
    # Script-specific command building
    if script_name == "Benchmark GREAC":
        # benchmark.sh: TRAIN TESTDIR GROUPNAME WINDOW METRIC KMER THRESHOLD REFERENCE
        cmd.extend([
            param_values.get("train_dir", ""),
            param_values.get("test_dir", ""),
            param_values.get("group_name", ""),
            param_values.get("window_size", ""),
            param_values.get("metric", ""),
            str(param_values.get("kmer_size", "")),
            str(param_values.get("threshold", "")),
            param_values.get("reference_path", "")
        ])
    
    elif script_name == "FASTA Regions":
        # fasta-regions.sh: GROUPNAME WINDOW [INPUT]
        cmd.extend([
            param_values.get("group_name", ""),
            param_values.get("window_size", "")
        ])
        dataset = param_values.get("dataset", "")
        if dataset and dataset.strip():
            cmd.append(dataset)
    
    elif script_name == "Feature Extraction":
        # extract-features.sh: TRAIN GROUPNAME WINDOW METRIC THRESHOLD REFERENCE
        cmd.extend([
            param_values.get("train_dir", ""),
            param_values.get("group_name", ""),
            param_values.get("window_size", ""),
            param_values.get("metric", ""),
            str(param_values.get("threshold", "")),
            param_values.get("reference_path", "")
        ])
    
    elif script_name == "Classification":
        # file-classification.sh: TRAIN TESTFILE GROUPNAME WINDOW METRIC KMER THRESHOLD REFERENCE
        cmd.extend([
            param_values.get("train_dir", ""),
            param_values.get("test_file", ""),
            param_values.get("group_name", ""),
            param_values.get("window_size", ""),
            param_values.get("metric", ""),
            str(param_values.get("kmer_size", "")),
            str(param_values.get("threshold", "")),
            param_values.get("reference_path", "")
        ])
    
    elif script_name == "Complete Pipeline (doall.sh)":
        # doall.sh: GROUPNAME WINDOW KMERSIZE THRESHOLD
        cmd.extend([
            param_values.get("group_name", ""),
            param_values.get("window_size", ""),
            str(param_values.get("kmer_size", "")),
            str(param_values.get("threshold", ""))
        ])
    
    return cmd

def expand_path_pattern(pattern: str) -> List[str]:
    """Expands patterns with wildcards and ~ to list of real paths"""
    # Expand ~ to home directory
    expanded_pattern = os.path.expanduser(pattern)
    
    # If there are no wildcards, return as list
    if '*' not in expanded_pattern and '?' not in expanded_pattern:
        return [expanded_pattern] if os.path.exists(expanded_pattern) else []
    
    # Use glob to expand wildcards
    matched_paths = glob(expanded_pattern, recursive=False)
    
    # Filter only directories or files that exist
    return [p for p in matched_paths if os.path.exists(p)]

def find_output_files() -> Dict[str, str]:
    """Searches for output files in configured directories/patterns"""
    found_items = {}
    
    for pattern in OUTPUT_DIRECTORIES:
        expanded_paths = expand_path_pattern(pattern)
        
        for path in expanded_paths:
            if os.path.isfile(path):
                # If it's a file, add directly
                file_info = f"📄 {os.path.basename(path)} ({format_file_size(os.path.getsize(path))})\n"
                file_info += f"Path: {path}"
                found_items[path] = file_info
            elif os.path.isdir(path):
                # If it's a directory, list the tree
                if os.listdir(path):  # Only add if not empty
                    found_items[path] = list_directory_tree(path)
    
    return found_items

def create_results_zip(output_paths: List[str]) -> str:
    """Creates a ZIP with all found results"""
    zip_path = "/tmp/all_results.zip"
    
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for path in output_paths:
            expanded_path = os.path.expanduser(path)
            if os.path.isfile(expanded_path):
                # Add individual file
                zipf.write(expanded_path, os.path.basename(expanded_path))
            elif os.path.isdir(expanded_path):
                # Add directory recursively
                for root, dirs, files in os.walk(expanded_path):
                    for file in files:
                        file_path = os.path.join(root, file)
                        arcname = os.path.join(
                            os.path.basename(expanded_path),
                            os.path.relpath(file_path, expanded_path)
                        )
                        zipf.write(file_path, arcname)
    
    return zip_path

# ============================================================================
# MAIN INTERFACE
# ============================================================================

def main():
    st.title("🧬 GREAC - Genomic Region Extraction and Classifier")
    st.markdown("---")
    
    # Script Selection
    st.header("📋 Script Selection")
    
    script_names = list(SCRIPTS_CONFIG.keys())
    selected_script = st.selectbox(
        "Choose the script to execute:",
        options=script_names,
        key="script_selector"
    )
    
    if selected_script:
        script_config = SCRIPTS_CONFIG[selected_script]
        
        # Show script description
        st.info(f"ℹ️ {script_config['description']}")
        
        # Check if script exists
        script_path = os.path.expanduser(script_config["script_path"])
        
        if not os.path.exists(script_path):
            st.warning(f"⚠️ Script not found at: {script_path}")
        
        # Render parameters dynamically
        st.header("⚙️ Script Parameters")
        st.markdown("*Fields marked with * are required*")
        
        param_values = {}
        parameters = script_config["parameters"]
        
        # Render parameters in two passes to allow dependencies
        # First pass: render independent parameters
        for param in parameters:
            if param.get("type") not in ["window_selector"]:
                param_values[param["name"]] = render_parameter(
                    param, 
                    f"{selected_script}_param"
                )
        
        # Second pass: render dependent parameters (window_selector needs group)
        for param in parameters:
            if param.get("type") == "window_selector":
                context = {"group_name": param_values.get("group_name")}
                param_values[param["name"]] = render_parameter(
                    param,
                    f"{selected_script}_param",
                    context=context
                )
        
        # Required field validation
        missing_required = []
        for param in parameters:
            if param.get("required"):
                value = param_values.get(param["name"])
                if not value or (isinstance(value, str) and not value.strip()):
                    missing_required.append(param["label"])
        
        st.markdown("---")
        
        # Execution Section
        st.header("🚀 Execution")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            if missing_required:
                st.error(f"⚠️ Missing required fields: {', '.join(missing_required)}")
                st.button("▶️ Start Processing", disabled=True)
            else:
                if st.button("▶️ Start Processing", type="primary", key="start_btn"):
                    # Build command
                    cmd = build_command(selected_script, script_path, parameters, param_values)
                    
                    st.info("📡 Starting execution...")
                    st.code(" ".join(cmd), language="bash")
                    
                    # Execute the script
                    with st.spinner("🔄 Processing..."):
                        try:
                            # Make script executable
                            os.chmod(script_path, 0o755)
                            
                            process = subprocess.Popen(
                                cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                universal_newlines=True,
                                shell=False
                            )
                            
                            output_placeholder = st.empty()
                            output_text = ""
                            
                            for line in process.stdout:
                                output_text += line
                                output_placeholder.code(output_text, language="bash")
                            
                            process.wait()
                            
                            if process.returncode == 0:
                                st.success("✅ Processing completed successfully!")
                            else:
                                st.error(f"❌ Process finished with error code: {process.returncode}")
                                
                        except Exception as e:
                            st.error(f"❌ Error during execution: {str(e)}")
        
        with col2:
            st.subheader("📥 Download Results")
            
            # Search for outputs
            found_outputs = find_output_files()
            
            if found_outputs:
                st.success(f"✅ Found {len(found_outputs)} director(y/ies) with results")
                
                if st.button("💾 Download All Results", key="download_all_btn"):
                    try:
                        zip_path = create_results_zip(list(found_outputs.keys()))
                        
                        with open(zip_path, "rb") as f:
                            st.download_button(
                                label="📦 Download ZIP",
                                data=f.read(),
                                file_name=f"{selected_script.replace(' ', '_')}_results.zip",
                                mime="application/zip",
                                key="download_zip_btn"
                            )
                    except Exception as e:
                        st.error(f"❌ Error creating ZIP: {str(e)}")
            else:
                st.warning("⚠️ No results found in configured directories")
        
        # Output Visualization Section
        st.markdown("---")
        st.header("📂 Output Files")
        
        found_outputs = find_output_files()
        
        if found_outputs:
            # Separate files and directories
            files_found = {k: v for k, v in found_outputs.items() if os.path.isfile(os.path.expanduser(k))}
            dirs_found = {k: v for k, v in found_outputs.items() if os.path.isdir(os.path.expanduser(k))}
            
            # Show statistics
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Items", len(found_outputs))
            with col2:
                st.metric("Directories", len(dirs_found))
            with col3:
                st.metric("Files", len(files_found))
            
            # Create tabs for files and directories
            if files_found and dirs_found:
                tab1, tab2 = st.tabs(["📁 Directories", "📄 Individual Files"])
                
                with tab1:
                    dir_tabs = st.tabs([os.path.basename(d) for d in dirs_found.keys()])
                    for idx, (output_dir, tree) in enumerate(dirs_found.items()):
                        with dir_tabs[idx]:
                            st.text(f"Path: {output_dir}")
                            st.code(tree, language="text")
                
                with tab2:
                    for file_path, file_info in files_found.items():
                        with st.expander(f"📄 {os.path.basename(file_path)}"):
                            st.code(file_info, language="text")
                            
                            # Individual download button
                            try:
                                expanded_path = os.path.expanduser(file_path)
                                with open(expanded_path, "rb") as f:
                                    st.download_button(
                                        label="⬇️ Download file",
                                        data=f.read(),
                                        file_name=os.path.basename(file_path),
                                        key=f"download_{file_path}"
                                    )
                            except Exception as e:
                                st.error(f"Error reading file: {str(e)}")
            
            elif dirs_found:
                dir_tabs = st.tabs([os.path.basename(d) for d in dirs_found.keys()])
                for idx, (output_dir, tree) in enumerate(dirs_found.items()):
                    with dir_tabs[idx]:
                        st.text(f"Path: {output_dir}")
                        st.code(tree, language="text")
            
            elif files_found:
                for file_path, file_info in files_found.items():
                    with st.expander(f"📄 {os.path.basename(file_path)}"):
                        st.code(file_info, language="text")
                        
                        # Individual download button
                        try:
                            expanded_path = os.path.expanduser(file_path)
                            with open(expanded_path, "rb") as f:
                                st.download_button(
                                    label="⬇️ Download file",
                                    data=f.read(),
                                    file_name=os.path.basename(file_path),
                                    key=f"download_{file_path}"
                                )
                        except Exception as e:
                            st.error(f"Error reading file: {str(e)}")
        else:
            st.info("ℹ️ No output files found. Run a script to generate results.")
            st.markdown("**Monitored patterns:**")
            for pattern in OUTPUT_DIRECTORIES:
                expanded = expand_path_pattern(pattern)
                if expanded:
                    st.success(f"✅ {pattern}")
                    for path in expanded:
                        st.text(f"    → {path}")
                else:
                    st.text(f"  • {pattern} (not found)")


if __name__ == "__main__":
    main()