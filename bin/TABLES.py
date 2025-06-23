import pandas as pd
import matplotlib.pyplot as plt
from textwrap import fill

# Carrega os dados
df = pd.read_excel('TABELAS.xlsx', index_col=0)

# Cores
cell_color = '#C7F464'
header_color = '#A5D200'

for col in df.columns:
    fig, ax = plt.subplots(figsize=(5, 6))
    ax.axis('off')

    # Cria a primeira linha como cabeçalho mesclado
    table_data = [[col, '']]  # Cabeçalho que simula colspan
    for row in df.index:
        label = fill(str(row), width=15)
        value = fill(str(df.loc[row, col]), width=25)
        table_data.append([label, value])

    # Cria tabela
    table = ax.table(
        cellText=table_data,
        cellLoc='left',
        loc='center',
        colWidths=[0.5, 0.5]
    )

    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2.5)

    for (row_idx, col_idx), cell in table.get_celld().items():
        if row_idx == 0:
            cell.set_facecolor(header_color)
            cell.set_text_props(ha='center', va='center', fontweight='bold')
            if col_idx == 1:
                cell.set_visible(False)  # Oculta a segunda célula da linha 0 (simulando colspan)
        else:
            cell.set_facecolor(cell_color)
            cell.set_text_props(va='center', ha='left')
            if col_idx == 0:
                cell.get_text().set_fontweight('bold')
 
   # Salvar imagem
    filename = f'tabela_{col}.png'.replace(" ", "_")
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Figura salva: {filename}")
    plt.close()
