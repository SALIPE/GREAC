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
# CONFIGURAÇÃO DOS SCRIPTS
# ============================================================================
GROUPS_BASE_PATH = "~/.project_cache"
SCRIPTS_CONFIG = {
    "Benchmark GREAC": {
        "script_path": "../scripts/local/benchmark.sh",
        "description": "Executa o benchmark completo do GREAC",
        "parameters": [
            {
                "name": "train_dir",
                "label": "📂 Diretório de Treino",
                "type": "text",
                "default": "~/Desktop/datasets/dengue/train/kmers",
                "required": True
            },
            {
                "name": "test_dir",
                "label": "📂 Diretório de Teste",
                "type": "text",
                "default": "~/Desktop/datasets/dengue/test",
                "required": True
            },
            {
                "name": "group_name",
                "label": "📁 Nome do Grupo",
                "type": "group_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "window_size",
                "label": "🪟 Janela/Window Size",
                "type": "window_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "metric",
                "label": "Métrica",
                "type": "select",
                "options": ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"],
                "default": "manhattan",
                "required": True
            },
            {
                "name": "kmer_size",
                "label": "Tamanho do K-mer",
                "type": "number",
                "default": 0,
                "step": 1,
                "required": True
            },
            {
                "name": "threshold",
                "label": "Percentual de consideração de janela",
                "type": "number",
                "default": 0.0,
                "step": 0.1,
                "format": "%.1f",
                "required": True
            },
             {
                "name": "reference_path",
                "label": "📄 Arquivo de Referência (.fasta)",
                "type": "file_path",
                "extensions": [".fasta", ".fa", ".fas"],
                "default": "~/Desktop/datasets/denv/",
                "required": True
            },
            {
                "name": "cache",
                "label": "Usar Cache",
                "type": "checkbox",
                "default": True,
                "required": False
            }
        ]
    },
    "Regiões FASTA": {
        "script_path": "../scripts/local/fasta-regions.sh",
        "description": "Retorna o .bed file com as regiões e reduz o dataset (se passado)",
        "parameters": [
            {
                "name": "group_name",
                "label": "📁 Nome do Grupo",
                "type": "group_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "window_size",
                "label": "🪟 Janela/Window Size",
                "type": "window_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "dataset",
                "label": "📂 Diretório da referencia",
                "type": "text",
                "default": None,
                "required": False
            }
        ]
    },
    "Extração de Features": {
        "script_path": "../scripts/local/extract-features.sh",
        "description": "Aplica o metodo do GREAC para extração de regiões e treinamento do modelo",
        "parameters": [
            {
                "name": "train_dir",
                "label": "📂 Diretório de Treino",
                "type": "text",
                "default": "~/Desktop/datasets/dengue/train/kmers",
                "required": True
            },
            {
                "name": "group_name",
                "label": "📁 Nome do Grupo",
                "type": "group_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "window_size",
                "label": "🪟 Janela/Window Size",
                "type": "window_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "metric",
                "label": "Métrica",
                "type": "select",
                "options": ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"],
                "default": "manhattan",
                "required": True
            },
            {
                "name": "threshold",
                "label": "Percentual de consideração de janela",
                "type": "number",
                "default": 0.0,
                "step": 0.1,
                "format": "%.1f",
                "required": True
            },
             {
                "name": "reference_path",
                "label": "📄 Arquivo de Referência (.fasta)",
                "type": "file_path",
                "extensions": [".fasta", ".fa", ".fas"],
                "default": "~/Desktop/datasets/denv/",
                "required": True
            },
            {
                "name": "cache",
                "label": "Usar Cache",
                "type": "checkbox",
                "default": True,
                "required": False
            }
        ]
    },
    "Classificação": {
        "script_path": "../scripts/local/file-classification.sh",
        "description": "Classifica um arquivo baseado em um modelo pré-treinado",
        "parameters": [
            {
                "name": "train_dir",
                "label": "📂 Diretório de Treino",
                "type": "text",
                "default": "~/Desktop/datasets/dengue/train/kmers",
                "required": True
            },
            {
                "name": "test_dir",
                "label": "📂 Diretório de Teste",
                "type": "text",
                "default": "~/Desktop/datasets/dengue/test",
                "required": True
            },
            {
                "name": "group_name",
                "label": "📁 Nome do Grupo",
                "type": "group_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "window_size",
                "label": "🪟 Janela/Window Size",
                "type": "window_selector",
                "base_path": "~/.project_cache",
                "default": "",
                "required": True
            },
            {
                "name": "metric",
                "label": "Métrica",
                "type": "select",
                "options": ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"],
                "default": "manhattan",
                "required": True
            },
            {
                "name": "kmer_size",
                "label": "Tamanho do K-mer",
                "type": "number",
                "default": 0,
                "step": 1,
                "required": True
            },
            {
                "name": "threshold",
                "label": "Percentual de consideração de janela",
                "type": "number",
                "default": 0.0,
                "step": 0.1,
                "format": "%.1f",
                "required": True
            },
             {
                "name": "reference_path",
                "label": "📄 Arquivo de Referência (.fasta)",
                "type": "file_path",
                "extensions": [".fasta", ".fa", ".fas"],
                "default": "~/Desktop/datasets/denv/",
                "required": True
            },
            {
                "name": "cache",
                "label": "Usar Cache",
                "type": "checkbox",
                "default": True,
                "required": False
            }
        ]
    }
}

# Diretórios/padrões onde os outputs podem ser encontrados (suporta wildcards)
OUTPUT_DIRECTORIES = [
    "~/.project_cache",
    "../GREAC/output-*",
    "../GREAC/*.pdf",
]

# ============================================================================
# FUNÇÕES AUXILIARES
# ============================================================================
def get_available_groups(base_path: str) -> List[str]:
    """Lista diretórios disponíveis como grupos"""
    if not os.path.exists(base_path):
        return []
    
    try:
        groups = [d for d in os.listdir(base_path) 
                 if os.path.isdir(os.path.join(base_path, d)) and not d.startswith('.')]
        return sorted(groups)
    except:
        return []

def get_available_windows(base_path: str, group_name: str) -> List[str]:
    """Lista subdiretórios disponíveis como janelas dentro de um grupo"""
    if not group_name or not os.path.exists(base_path):
        return []
    
    group_path = os.path.join(base_path, group_name)
    if not os.path.exists(group_path):
        return []
    
    try:
        windows = [d for d in os.listdir(group_path) 
                  if os.path.isdir(os.path.join(group_path, d)) and not d.startswith('.')]
        return sorted(windows)
    except:
        return []
    
def list_directory_tree(directory: str, prefix: str = "", max_depth: int = 5) -> str:
    """Lista a estrutura de diretórios em formato de árvore"""
    if not os.path.exists(directory):
        return f"❌ Diretório não encontrado: {directory}\n"
    
    tree_str = f"📁 {os.path.basename(directory)}/\n"
    
    try:
        items = []
        for root, dirs, files in os.walk(directory):
            level = root.replace(directory, "").count(os.sep)
            if level >= max_depth:
                continue
                
            indent = "    " * level
            
            # Adiciona diretórios
            for d in sorted(dirs):
                items.append(f"{indent}├── 📁 {d}/")
            
            # Adiciona arquivos
            for f in sorted(files):
                file_path = os.path.join(root, f)
                size = os.path.getsize(file_path)
                size_str = format_file_size(size)
                items.append(f"{indent}├── 📄 {f} ({size_str})")
        
        tree_str += "\n".join(items)
    except Exception as e:
        tree_str += f"❌ Erro ao listar: {str(e)}"
    
    return tree_str

def format_file_size(size_bytes: int) -> str:
    """Formata o tamanho do arquivo em formato legível"""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"

def render_parameter(param_config: Dict[str, Any], key_prefix: str, context: Dict[str, Any] = None) -> Any:
    """Renderiza um parâmetro de acordo com seu tipo"""
    param_type = param_config.get("type", "text")
    name = param_config["name"]
    label = param_config.get("label", name)
    default = param_config.get("default", "")
    required = param_config.get("required", False)
    
    key = f"{key_prefix}_{name}"
    
    if param_type == "text":
        value = st.text_input(
            label + (" *" if required else ""),
            value=default,
            key=key
        )
    elif param_type == "group_selector":
        # Campo híbrido para seleção/digitação de grupo
        base_path = os.path.expanduser(GROUPS_BASE_PATH)
        available_groups = get_available_groups(base_path)
        
        st.markdown(f"**{label}**" + (" *" if required else ""))
        
        col1, col2 = st.columns([3, 1])
        with col1:
            if available_groups:
                # Radio button para escolher modo
                mode = st.radio(
                    "Modo de entrada:",
                    options=["Selecionar existente", "Digitar novo"],
                    key=f"{key}_mode",
                    horizontal=True
                )
                
                if mode == "Selecionar existente":
                    value = st.selectbox(
                        "Escolha um grupo:",
                        options=[""] + available_groups,
                        key=f"{key}_select"
                    )
                else:
                    value = st.text_input(
                        "Digite o nome do grupo:",
                        value=default,
                        key=f"{key}_text"
                    )
            else:
                st.info(f"ℹ️ Nenhum grupo encontrado em: {base_path}")
                value = st.text_input(
                    "Digite o nome do grupo:",
                    value=default,
                    key=f"{key}_text"
                )
        
        with col2:
            if available_groups:
                st.metric("Grupos disponíveis", len(available_groups))
    
    elif param_type == "window_selector":
        # Campo híbrido para seleção/digitação de janela (depende do grupo)
        base_path =  os.path.expanduser(GROUPS_BASE_PATH)
        group_name = context.get("group_name") if context else None
        
        st.markdown(f"**{label}**" + (" *" if required else ""))
        
        if not group_name:
            st.warning("⚠️ Selecione um grupo primeiro para ver janelas disponíveis")
            value = st.text_input(
                "Digite o valor:",
                value=default,
                key=f"{key}_text"
            )
        else:
            available_windows = get_available_windows(base_path, group_name)
            
            col1, col2 = st.columns([3, 1])
            with col1:
                if available_windows:
                    # Radio button para escolher modo
                    mode = st.radio(
                        "Modo de entrada:",
                        options=["Selecionar existente", "Digitar novo"],
                        key=f"{key}_mode",
                        horizontal=True
                    )
                    
                    if mode == "Selecionar existente":
                        value = st.selectbox(
                            f"Escolha uma janela (do grupo '{group_name}'):",
                            options=[""] + available_windows,
                            key=f"{key}_select"
                        )
                    else:
                        value = st.text_input(
                            "Digite o valor:",
                            value=default,
                            key=f"{key}_text"
                        )
                else:
                    st.info(f"ℹ️ Nenhuma janela encontrada para o grupo '{group_name}'")
                    value = st.text_input(
                        "Digite o valor:",
                        value=default,
                        key=f"{key}_text"
                    )
            
            with col2:
                if available_windows:
                    st.metric("Janelas disponíveis", len(available_windows))
    
    elif param_type == "number":
        step = param_config.get("step", 1)
        format_str = param_config.get("format", "%d")
        value = st.number_input(
            label + (" *" if required else ""),
            value=default,
            step=step,
            format=format_str,
            key=key
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
            key=key
        )
    elif param_type == "checkbox":
        value = st.checkbox(
            label,
            value=default,
            key=key
        )
    elif param_type == "file_path":
        value = st.text_input(
            label + (" *" if required else ""),
            value=default,
            key=key,
            help=f"Extensões aceitas: {', '.join(param_config.get('extensions', []))}"
        )
    else:
        value = st.text_input(label, value=default, key=key)
    
    return value

def build_command(script_path: str, parameters: List[Dict], param_values: Dict) -> List[str]:
    """Constrói o comando a ser executado"""
    cmd = [script_path]
    
    for param in parameters:
        name = param["name"]
        value = param_values.get(name)
        
        if param.get("type") == "checkbox":
            if not value and name == "cache":
                cmd.append("--no-cache")
        else:
            if value is not None and str(value).strip():
                cmd.append(str(value))
    
    return cmd

def expand_path_pattern(pattern: str) -> List[str]:
    """Expande padrões com wildcards e ~ para lista de caminhos reais"""
    # Expande ~ para home directory
    expanded_pattern = os.path.expanduser(pattern)
    
    # Se não tem wildcards, retorna como lista
    if '*' not in expanded_pattern and '?' not in expanded_pattern:
        return [expanded_pattern] if os.path.exists(expanded_pattern) else []
    
    # Usa glob para expandir wildcards
    matched_paths = glob(expanded_pattern, recursive=False)
    
    # Filtra apenas diretórios ou arquivos que existem
    return [p for p in matched_paths if os.path.exists(p)]

def find_output_files() -> Dict[str, str]:
    """Procura por arquivos de output nos diretórios/padrões configurados"""
    found_items = {}
    
    for pattern in OUTPUT_DIRECTORIES:
        expanded_paths = expand_path_pattern(pattern)
        
        for path in expanded_paths:
            if os.path.isfile(path):
                # Se é um arquivo, adiciona diretamente
                file_info = f"📄 {os.path.basename(path)} ({format_file_size(os.path.getsize(path))})\n"
                file_info += f"Caminho: {path}"
                found_items[path] = file_info
            elif os.path.isdir(path):
                # Se é diretório, lista a árvore
                if os.listdir(path):  # Só adiciona se não estiver vazio
                    found_items[path] = list_directory_tree(path)
    
    return found_items

def create_results_zip(output_paths: List[str]) -> str:
    """Cria um ZIP com todos os resultados encontrados"""
    zip_path = "/tmp/all_results.zip"
    
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for path in output_paths:
            if os.path.isfile(path):
                # Adiciona arquivo individual
                zipf.write(path, os.path.basename(path))
            elif os.path.isdir(path):
                # Adiciona diretório recursivamente
                for root, dirs, files in os.walk(path):
                    for file in files:
                        file_path = os.path.join(root, file)
                        arcname = os.path.join(
                            os.path.basename(path),
                            os.path.relpath(file_path, path)
                        )
                        zipf.write(file_path, arcname)
    
    return zip_path

# ============================================================================
# INTERFACE PRINCIPAL
# ============================================================================

def main():
    st.title("🧬 GREAC - Genomic Region Extraction and Classifier")
    st.markdown("---")
    
    # Seleção do Script
    st.header("📋 Seleção de Script")
    
    script_names = list(SCRIPTS_CONFIG.keys())
    selected_script = st.selectbox(
        "Escolha o script a executar:",
        options=script_names,
        key="script_selector"
    )
    
    if selected_script:
        script_config = SCRIPTS_CONFIG[selected_script]
        
        # Mostra descrição do script
        st.info(f"ℹ️ {script_config['description']}")
        
        # Verifica se o script existe
        script_path = script_config["script_path"]
        
       # Renderiza parâmetros dinamicamente
        st.header("⚙️ Parâmetros do Script")
        st.markdown("*Campos marcados com * são obrigatórios*")
        
        param_values = {}
        parameters = script_config["parameters"]
        
        # Renderiza parâmetros em duas passadas para permitir dependências
        # Primeira passada: renderiza parâmetros independentes
        for param in parameters:
            if param.get("type") not in ["window_selector"]:
                param_values[param["name"]] = render_parameter(
                    param, 
                    f"{selected_script}_param"
                )
        
        # Segunda passada: renderiza parâmetros dependentes (window_selector precisa do grupo)
        for param in parameters:
            if param.get("type") == "window_selector":
                context = {"group_name": param_values.get("group_name")}
                param_values[param["name"]] = render_parameter(
                    param,
                    f"{selected_script}_param",
                    context=context
                )
        
        # Validação de campos obrigatórios
        missing_required = []
        for param in parameters:
            if param.get("required"):
                value = param_values.get(param["name"])
                if not value or (isinstance(value, str) and not value.strip()):
                    missing_required.append(param["label"])
        
        st.markdown("---")
        
        # Seção de Execução
        st.header("🚀 Execução")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            if missing_required:
                st.error(f"⚠️ Campos obrigatórios faltando: {', '.join(missing_required)}")
                st.button("▶️ Iniciar Processamento", disabled=True)
            else:
                if st.button("▶️ Iniciar Processamento", type="primary", key="start_btn"):
                    # Constrói comando
                    cmd = build_command(script_path, parameters, param_values)
                    
                    st.info("📡 Iniciando execução...")
                    st.code(" ".join(cmd), language="bash")
                    
                    # Executa o script
                    with st.spinner("🔄 Processando..."):
                        try:
                            process = subprocess.Popen(
                                cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                universal_newlines=True
                            )
                            
                            output_placeholder = st.empty()
                            output_text = ""
                            
                            for line in process.stdout:
                                output_text += line
                                output_placeholder.code(output_text, language="bash")
                            
                            process.wait()
                            
                            if process.returncode == 0:
                                st.success("✅ Processamento concluído com sucesso!")
                            else:
                                st.error(f"❌ Processo finalizado com código de erro: {process.returncode}")
                                
                        except Exception as e:
                            st.error(f"❌ Erro durante a execução: {str(e)}")
        
        with col2:
            st.subheader("📥 Download de Resultados")
            
            # Procura por outputs
            found_outputs = find_output_files()
            
            if found_outputs:
                st.success(f"✅ Encontrados {len(found_outputs)} diretório(s) com resultados")
                
                if st.button("💾 Baixar Todos os Resultados", key="download_all_btn"):
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
                        st.error(f"❌ Erro ao criar ZIP: {str(e)}")
            else:
                st.warning("⚠️ Nenhum resultado encontrado nos diretórios configurados")
        
        # Seção de Visualização de Outputs
        st.markdown("---")
        st.header("📂 Arquivos de Output")
        
        found_outputs = find_output_files()
        
        if found_outputs:
            # Separa arquivos e diretórios
            files_found = {k: v for k, v in found_outputs.items() if os.path.isfile(k)}
            dirs_found = {k: v for k, v in found_outputs.items() if os.path.isdir(k)}
            
            # Mostra estatísticas
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total de Itens", len(found_outputs))
            with col2:
                st.metric("Diretórios", len(dirs_found))
            with col3:
                st.metric("Arquivos", len(files_found))
            
            # Cria tabs para arquivos e diretórios
            if files_found and dirs_found:
                tab1, tab2 = st.tabs(["📁 Diretórios", "📄 Arquivos Individuais"])
                
                with tab1:
                    dir_tabs = st.tabs([os.path.basename(d) for d in dirs_found.keys()])
                    for idx, (output_dir, tree) in enumerate(dirs_found.items()):
                        with dir_tabs[idx]:
                            st.text(f"Caminho: {output_dir}")
                            st.code(tree, language="text")
                
                with tab2:
                    for file_path, file_info in files_found.items():
                        with st.expander(f"📄 {os.path.basename(file_path)}"):
                            st.code(file_info, language="text")
                            
                            # Botão de download individual
                            try:
                                with open(file_path, "rb") as f:
                                    st.download_button(
                                        label="⬇️ Baixar arquivo",
                                        data=f.read(),
                                        file_name=os.path.basename(file_path),
                                        key=f"download_{file_path}"
                                    )
                            except Exception as e:
                                st.error(f"Erro ao ler arquivo: {str(e)}")
            
            elif dirs_found:
                dir_tabs = st.tabs([os.path.basename(d) for d in dirs_found.keys()])
                for idx, (output_dir, tree) in enumerate(dirs_found.items()):
                    with dir_tabs[idx]:
                        st.text(f"Caminho: {output_dir}")
                        st.code(tree, language="text")
            
            elif files_found:
                for file_path, file_info in files_found.items():
                    with st.expander(f"📄 {os.path.basename(file_path)}"):
                        st.code(file_info, language="text")
                        
                        # Botão de download individual
                        try:
                            with open(file_path, "rb") as f:
                                st.download_button(
                                    label="⬇️ Baixar arquivo",
                                    data=f.read(),
                                    file_name=os.path.basename(file_path),
                                    key=f"download_{file_path}"
                                )
                        except Exception as e:
                            st.error(f"Erro ao ler arquivo: {str(e)}")
        else:
            st.info("ℹ️ Nenhum arquivo de output encontrado. Execute um script para gerar resultados.")
            st.markdown("**Padrões monitorados:**")
            for pattern in OUTPUT_DIRECTORIES:
                expanded = expand_path_pattern(pattern)
                if expanded:
                    st.success(f"✅ {pattern}")
                    for path in expanded:
                        st.text(f"    → {path}")
                else:
                    st.text(f"  • {pattern} (não encontrado)")


if __name__ == "__main__":
    main()