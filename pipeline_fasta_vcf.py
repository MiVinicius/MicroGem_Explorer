import os
import subprocess
import sys

# Configs
ref = "ref.fasta"
sample = "sample.fasta"
vcf_out = "resultado.vcf.gz"
env_name = "fasta2vcf"

# Conda
def check_conda():
    try:
        subprocess.run(["conda", "--version"], check=True)
        print("Conda/Miniforge detectado")
    except:
        print("Conda/Miniforge não encontrado. Baixe e instale:")
        print("https://github.com/conda-forge/miniforge/releases")
        sys.exit(1)

# Criar ambiente e instalar ferramentas
def criar_ambiente():
    print(f"Criando ambiente Conda '{env_name}' com minimap2, samtools e bcftools...")
    subprocess.run([
        "conda", "create", "-y", "-n", env_name,
        "minimap2", "samtools", "bcftools", "python"
    ], check=True)
    print("Ambiente criado e ferramentas instaladas")

# Função para rodar comando dentro do ambiente Conda
def run_in_env(cmd):
    #full_cmd = f"conda run -n {env_name} " + cmd
    subprocess.run(cmd, shell=True, check=True)

# Pipeline FASTA → BAM → VCF
def fasta_para_vcf():
    print("Indexando referência...")
    run_in_env(f"samtools faidx {ref}")

    print("Alinhando amostra contra referência...")
    run_in_env(f"minimap2 -ax asm5 {ref} {sample} > alinhamento.sam")

    print("Convertendo SAM → BAM...")
    run_in_env("samtools view -S -b alinhamento.sam > alinhamento.bam")

    print("Ordenando BAM...")
    run_in_env("samtools sort alinhamento.bam -o alinhamento.sorted.bam")

    print("Indexando BAM...")
    run_in_env("samtools index alinhamento.sorted.bam")

    print("Chamando variantes...")
    run_in_env(f"bcftools mpileup -f {ref} alinhamento.sorted.bam | bcftools call -mv -Oz -o {vcf_out}")

    print("Indexando VCF...")
    run_in_env(f"bcftools index {vcf_out}")

    print(f"Pipeline concluído! VCF gerado: {vcf_out}")

# Execução principal
if __name__ == "__main__":
    # Checa arquivos
    if not os.path.exists(ref) or not os.path.exists(sample):
        print("Verifique se os arquivos FASTA existem no diretório atual")
        sys.exit(1)

    check_conda()
    #criar_ambiente()
    fasta_para_vcf()
