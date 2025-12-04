import pandas as pd
import csv

class GFFParser:
    def __init__(self, gff_path):
        self.gff_path = gff_path
        self.genes = []

    def parse(self):
        """Analisa o arquivo GFF e extrai características de genes."""
        try:
            with open(self.gff_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\t')
                for line in reader:
                    if not line or line[0].startswith('#'):
                        continue
                    
                    if len(line) < 9:
                        continue

                    # Colunas GFF3: seqid, source, type, start, end, score, strand, phase, attributes
                    feature_type = line[2]
                    
                    # Estamos interessados em características 'gene' ou 'CDS'
                    if feature_type not in ['gene', 'CDS']:
                        continue

                    chrom = line[0]
                    start = int(line[3])
                    end = int(line[4])
                    attributes_str = line[8]
                    
                    # Analisar atributos para obter ID ou Nome
                    attributes = {}
                    for attr in attributes_str.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attributes[key] = value
                    
                    gene_name = attributes.get('Name', attributes.get('ID', 'Unknown'))

                    self.genes.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'name': gene_name,
                        'type': feature_type
                    })
            print(f"Sucesso ao analisar {len(self.genes)} características de genes.")
        except Exception as e:
            print(f"Erro ao analisar GFF: {e}")
            raise

    def get_genes(self, chrom, pos):
        """Retorna uma lista de genes que se sobrepõem à posição dada."""
        overlapping_genes = []
        for gene in self.genes:
            if gene['chrom'] == chrom and gene['start'] <= pos <= gene['end']:
                overlapping_genes.append(gene['name'])
        return overlapping_genes
