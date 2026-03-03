import streamlit as st
import subprocess
import os
import shutil
import pandas as pd
from io import BytesIO
from Bio import Entrez
import plotly.express as px

# --- Настройки путей ---
DEFAULT_DIAMOND_PATH = r"C:\Users\diamond\diamond.exe"
TEMP_DIR = "diamond_antagonism_temp"

DIAMOND_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

def setup_temp_dir():
    if os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)
    os.makedirs(TEMP_DIR, exist_ok=True)

def run_command(command):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Ошибка выполнения: {' '.join(command)}\n{result.stderr}")
    return result.stdout

def build_database(diamond_path, reference_fasta_path, db_path):
    cmd = [diamond_path, "makedb", "--in", reference_fasta_path, "-d", db_path, "--quiet"]
    run_command(cmd)

def run_diamond_blastx(diamond_path, db_path, query_dna_path, output_tsv, threads, evalue):
    cmd = [
        diamond_path, "blastx", "-d", db_path, "-q", query_dna_path,
        "-o", output_tsv, "-f", "6", "--threads", str(threads),
        "--evalue", str(evalue), "--quiet"
    ]
    run_command(cmd)

def to_excel(df):
    """Конвертирует DataFrame в байты Excel-файла"""
    output = BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Results')
    processed_data = output.getvalue()
    return processed_data

# --- Интерфейс Streamlit ---
st.set_page_config(page_title="Genome Miner Pro", layout="wide")
st.title("🦠 Genome Miner Pro: Факторы антагонизма")

# Создаем две вкладки
tab_analysis, tab_ncbi = st.tabs(["🧬 Анализ геномов (DIAMOND)", "🌐 Парсер NCBI (Создание базы)"])

# ==========================================
# Вкладка 1: АНАЛИЗ ГЕНОМОВ
# ==========================================
with tab_analysis:
    st.markdown("Здесь вы можете сравнить ваши сборки геномов с базой известных генов.")
    
    st.sidebar.header("Настройки DIAMOND")
    diamond_exe = st.sidebar.text_input("Путь к diamond.exe", DEFAULT_DIAMOND_PATH)
    evalue_cutoff = st.sidebar.number_input("Порог E-value", value=1e-5, format="%e")
    pident_cutoff = st.sidebar.slider("Минимальная идентичность (%)", 0, 100, 30)
    threads = st.sidebar.number_input("Количество потоков", min_value=1, value=4)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("1. База маркеров (Белки)")
        st.info("Загрузите .faa файл с референсными белками (его можно создать во вкладке 'Парсер NCBI').")
        ref_file = st.file_uploader("Загрузить референсы", type=['fasta', 'faa', 'fa'], key="ref")

    with col2:
        st.subheader("2. Исследуемые геномы (ДНК)")
        st.info("Загрузите файлы сборок (3_contigs.fasta и т.д.). Можно выбрать несколько.")
        genome_files = st.file_uploader("Загрузить геномы", type=['fasta', 'fna', 'fa'], accept_multiple_files=True, key="genomes")

    if st.button("🚀 Запустить анализ", use_container_width=True, type="primary"):
        if not os.path.exists(diamond_exe):
            st.error(f"DIAMOND не найден по пути: {diamond_exe}")
        elif not ref_file or not genome_files:
            st.warning("Загрузите и базу маркеров, и исследуемые геномы.")
        else:
            try:
                with st.spinner("Создание базы данных DIAMOND..."):
                    setup_temp_dir()
                    ref_path = os.path.join(TEMP_DIR, "references.faa")
                    with open(ref_path, "wb") as f:
                        f.write(ref_file.getbuffer())
                    
                    db_path = os.path.join(TEMP_DIR, "ref_db")
                    build_database(diamond_exe, ref_path, db_path)
                    
                all_results = []
                progress_bar = st.progress(0)
                
                for i, genome in enumerate(genome_files):
                    with st.spinner(f"Анализ: {genome.name}..."):
                        genome_path = os.path.join(TEMP_DIR, f"query_{i}.fasta")
                        with open(genome_path, "wb") as f:
                            f.write(genome.getbuffer())
                        
                        output_tsv = os.path.join(TEMP_DIR, f"result_{i}.tsv")
                        run_diamond_blastx(diamond_exe, db_path, genome_path, output_tsv, threads, evalue_cutoff)
                        
                        if os.path.getsize(output_tsv) > 0:
                            df = pd.read_csv(output_tsv, sep="\t", names=DIAMOND_COLUMNS)
                            df.insert(0, "Genome", genome.name)
                            df = df[df['pident'] >= pident_cutoff]
                            all_results.append(df)
                            
                    progress_bar.progress((i + 1) / len(genome_files))
                
                if all_results:
                    final_df = pd.concat(all_results, ignore_index=True)
                    final_df.rename(columns={
                        "qseqid": "Контиг (ДНК)", "sseqid": "Найденный маркер",
                        "pident": "Идентичность (%)", "evalue": "E-value",
                        "length": "Длина совпадения", "bitscore": "Bit-score"
                    }, inplace=True)
                    
                    st.success("✅ Анализ успешно завершен!")
                    
                    # --- ГРАФИКИ ---
                    st.subheader("📊 Визуализация результатов")
                    col_g1, col_g2 = st.columns(2)
                    
                    with col_g1:
                        # График 1: Количество найденных генов по штаммам
                        hits_count = final_df['Genome'].value_counts().reset_index()
                        hits_count.columns = ['Штамм', 'Количество маркеров']
                        fig1 = px.bar(hits_count, x='Штамм', y='Количество маркеров', 
                                      title="Количество найденных маркеров по штаммам",
                                      color='Штамм')
                        st.plotly_chart(fig1, use_container_width=True)
                        
                    with col_g2:
                        # График 2: Идентичность vs Bit-score
                        fig2 = px.scatter(final_df, x='Идентичность (%)', y='Bit-score', 
                                          color='Genome', hover_data=['Найденный маркер', 'Контиг (ДНК)'],
                                          title="Качество выравнивания (Идентичность vs Score)")
                        st.plotly_chart(fig2, use_container_width=True)

                    # --- ТАБЛИЦА И СКАЧИВАНИЕ ---
                    st.subheader("📋 Сводная таблица")
                    st.dataframe(final_df[["Genome", "Контиг (ДНК)", "Найденный маркер", "Идентичность (%)", "E-value", "Длина совпадения"]])
                    
                    col_d1, col_d2 = st.columns(2)
                    with col_d1:
                        # Скачать CSV
                        csv = final_df.to_csv(index=False).encode('utf-8')
                        st.download_button("📥 Скачать CSV", data=csv, file_name='results.csv', mime='text/csv', use_container_width=True)
                    with col_d2:
                        # Скачать Excel
                        excel_data = to_excel(final_df)
                        st.download_button("📥 Скачать Excel (.xlsx)", data=excel_data, file_name='results.xlsx', mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet', use_container_width=True)
                else:
                    st.warning("Ничего не найдено. Попробуйте снизить 'Минимальную идентичность' или увеличить E-value.")
                    
            except Exception as e:
                st.error(f"Ошибка выполнения: {e}")


# ==========================================
# Вкладка 2: ПАРСЕР NCBI
# ==========================================
with tab_ncbi:
    st.markdown("""
    Здесь вы можете найти и скачать белки из базы **NCBI Protein**, чтобы создать свой файл референсов (базу маркеров).
    *Пример запроса:* `nuclease Bacillus[Organism]` или `dispersin B biofilm`.
    """)
    
    ncbi_email = st.text_input("Ваш Email (Обязательное требование NCBI)", placeholder="example@mail.com")
    search_term = st.text_input("Поисковый запрос", placeholder="Например: bacteriocin Pseudomonas")
    max_results = st.number_input("Максимум результатов для скачивания", min_value=1, max_value=500, value=20)
    
    if st.button("🔍 Искать и скачать последовательности"):
        if not ncbi_email:
            st.error("Пожалуйста, введите ваш Email. Это правило серверов NCBI для использования их API.")
        elif not search_term:
            st.warning("Введите поисковый запрос.")
        else:
            try:
                Entrez.email = ncbi_email
                with st.spinner("Поиск идентификаторов в NCBI..."):
                    handle = Entrez.esearch(db="protein", term=search_term, retmax=max_results)
                    record = Entrez.read(handle)
                    handle.close()
                    
                    id_list = record["IdList"]
                    count = record["Count"]
                    
                if not id_list:
                    st.warning(f"По запросу '{search_term}' ничего не найдено. Попробуйте изменить запрос.")
                else:
                    st.success(f"Найдено результатов всего: {count}. Скачиваем первые {len(id_list)}...")
                    
                    with st.spinner("Скачивание FASTA последовательностей..."):
                        fetch_handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text")
                        fasta_data = fetch_handle.read()
                        fetch_handle.close()
                        
                    st.subheader("Предпросмотр (первые 1000 символов):")
                    st.text(fasta_data[:1000] + "\n\n... (продолжение скрыто)")
                    
                    st.download_button(
                        label="💾 Сохранить базу как .faa (FASTA)",
                        data=fasta_data,
                        file_name=f"ncbi_database_{search_term.replace(' ', '_')}.faa",
                        mime="text/plain",
                        use_container_width=True,
                        type="primary"
                    )
                    st.info("👆 Скачайте этот файл и загрузите его во вкладке 'Анализ геномов' в качестве 'Базы маркеров'.")
                    
            except Exception as e:
                st.error(f"Произошла ошибка при обращении к NCBI: {e}")