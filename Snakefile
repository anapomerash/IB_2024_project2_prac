rule all:
    input:
        expand("results/csv/{sample}_parsed_variants.csv", sample=["SRR1705858", "SRR1705859", "SRR1705860"]),
        expand("results/vcf/{sample}_rare_variants_annotated.vcf", sample=["SRR1705858", "SRR1705859", "SRR1705860"]),
        "results/csv/significant_variants.csv",
        "results/logs/significant_variants_stats.txt",
        "results/vcf/roommate_annotated_variants.vcf",
        "results/logs/coverage_summary.txt",
        "results/csv/epitope_mutations.csv"

rule bwa_index:
    input:
        "data/reference.fasta"
    output:
        "data/reference.fasta.amb",
        "data/reference.fasta.ann",
        "data/reference.fasta.bwt",
        "data/reference.fasta.pac",
        "data/reference.fasta.sa"
    log:
        "results/logs/bwa_index.log"
    shell:
        """
        bwa index {input} > {log} 2>&1
        """

rule index_reference:
    input:
        "data/reference.fasta"
    output:
        "data/reference.fasta.fai"
    log:
        "results/logs/index_reference.log"
    shell:
        """
        samtools faidx {input} > {log} 2>&1
        """

rule calculate_read_statistics:
    input:
        "results/bam/roommate_sorted.bam"
    output:
        unmapped="results/logs/unmapped_reads.txt",
        total="results/logs/total_reads.txt"
    log:
        "results/logs/read_statistics.log"
    threads: 1
    shell:
        """
        samtools view -c -f4 {input} > {output.unmapped} 2>> {log}
        samtools view -c {input} > {output.total} 2>> {log}
        """

rule index_bam:
    input:
        "results/bam/roommate_sorted.bam"
    output:
        "results/bam/roommate_sorted.bam.bai"
    log:
        "results/logs/index_bam.log"
    shell:
        """
        samtools index {input} > {log} 2>&1
        """

rule calculate_depth_limit:
    input:
        ref_fai="data/reference.fasta.fai",
        total_reads="results/logs/total_reads.txt",
        fastq="data/roommate.fastq"
    output:
        "results/logs/depth_limit.txt",
        "results/logs/average_read_length.txt"
    log:
        "results/logs/calculate_depth_limit.log"
    threads: 1
    run:
        with open(input.ref_fai) as ref_fai:
            ref_length = int(ref_fai.readline().split()[1])

        with open(input.total_reads) as total_reads_file:
            total_reads = int(total_reads_file.read().strip())

        total_base_pairs = 0
        read_count = 0

        with open(input.fastq) as fastq_file:
            for i, line in enumerate(fastq_file):
                if i % 4 == 1:
                    total_base_pairs += len(line.strip())
                    read_count += 1

        average_read_length = total_base_pairs / read_count

        with open(output[1], "w") as avg_length_file:
            avg_length_file.write(f"{average_read_length:.2f}\n")

        coverage = (total_reads * average_read_length) / ref_length

        depth_limit = int(coverage * 1.2)

        with open(output[0], "w") as depth_limit_file:
            depth_limit_file.write(f"{depth_limit}\n")

rule create_mpileup:
    input:
        bam="results/bam/roommate_sorted.bam",
        bai="results/bam/roommate_sorted.bam.bai",
        ref="data/reference.fasta",
        depth_limit="results/logs/depth_limit.txt"
    output:
        "results/bam/roommate.mpileup"
    log:
        "results/logs/mpileup.log"
    threads: 4
    shell:
        """
        samtools mpileup -d $(cat {input.depth_limit}) -f {input.ref} {input.bam} > {output} 2> {log}
        """

rule run_varscan:
    input:
        mpileup="results/bam/roommate.mpileup"
    output:
        variants="results/vcf/roommate_high_freq_variants.vcf"
    log:
        "results/logs/varscan.log"
    shell:
        """
        java -jar /Users/anapomerash/bioinformatics_project/tools/VarScan.v2.3.9.jar mpileup2snp {input.mpileup} \
        --min-var-freq 0.95 \
        --output-vcf 1 > {output.variants} 2> {log}
        """

rule extract_variants:
    input:
        vcf="results/vcf/roommate_high_freq_variants.vcf"
    output:
        parsed_variants="results/vcf/roommate_parsed_variants.txt"
    log:
        "results/logs/extract_variants.log"
    shell:
        """
        awk 'NR>24 {{print $1, $2, $4, $5}}' {input.vcf} > {output.parsed_variants} 2> {log}
        """

rule annotate_variants:
    input:
        vcf="results/vcf/roommate_high_freq_variants.vcf"
    output:
        annotated_vcf="results/vcf/roommate_annotated_variants.vcf"
    log:
        "results/logs/annotate_variants.log"
    shell:
        """
        java -jar /Users/anapomerash/miniforge3/envs/myenv_x86/share/snpeff-5.2-1/snpEff.jar Flu {input.vcf} > {output.annotated_vcf} 2> {log}
        """

rule vcf_to_table:
    input:
        "results/vcf/roommate_annotated_variants.vcf"
    output:
        "results/vcf/roommate_annotated_table.csv"
    log:
        "results/logs/vcf_to_table.log"
    run:
        import pandas as pd
        import vcfpy

        # Чтение VCF файла
        vcf_file = input[0]
        output_file = output[0]

        # Парсинг VCF
        reader = vcfpy.Reader.from_path(vcf_file)
        data = []

        for record in reader:
            ann = record.INFO.get("ANN", ["N/A"])[0]
            data.append({
                "Chromosome": record.CHROM,
                "Position": record.POS,
                "ID": record.ID,
                "Reference": record.REF,
                "Alternate": ",".join([alt.value for alt in record.ALT]),
                "Quality": record.QUAL,
                "Filter": record.FILTER,
                "Annotation": ann
            })

        # Создание таблицы
        df = pd.DataFrame(data)

        # Сохранение в CSV
        df.to_csv(output_file, index=False)

        print(f"Таблица успешно создана: {output_file}")

rule run_varscan_rare:
    input:
        mpileup="results/bam/roommate.mpileup"
    output:
        variants="results/vcf/roommate_rare_variants.vcf"
    log:
        "results/logs/varscan_rare.log"
    shell:
        """
        java -jar /Users/anapomerash/bioinformatics_project/tools/VarScan.v2.3.9.jar mpileup2snp {input.mpileup} \
        --min-var-freq 0.001 \
        --output-vcf 1 > {output.variants} 2> {log}
        """

rule calculate_read_coverage:
    input:
        fastq=["data/SRR1705858.fastq", "data/SRR1705859.fastq", "data/SRR1705860.fastq"],
        reference="data/reference.fasta"
    output:
        "results/logs/coverage_summary.txt"
    log:
        "results/logs/calculate_coverage.log"
    run:
        import math

        # Функция для подсчёта ридов и средней длины
        def calculate_coverage(fastq_file, reference_length):
            total_reads = 0
            total_bases = 0

            with open(fastq_file, "r") as f:
                for i, line in enumerate(f):
                    if i % 4 == 1:  # Строка с последовательностью
                        total_reads += 1
                        total_bases += len(line.strip())

            avg_read_length = total_bases / total_reads if total_reads > 0 else 0
            coverage = (total_reads * avg_read_length) / reference_length
            return total_reads, avg_read_length, coverage

        # Получаем длину референса
        with open(f"{input.reference}.fai") as ref_fai:
            reference_length = int(ref_fai.readline().split()[1])

        # Вычисляем покрытие для каждого FASTQ
        with open(output[0], "w") as out, open(log[0], "w") as log_file:
            out.write(f"Reference length: {reference_length} bp\n")
            for fastq in input.fastq:
                total_reads, avg_read_length, coverage = calculate_coverage(fastq, reference_length)
                out.write(f"{fastq}: {total_reads} reads, average read length {avg_read_length:.2f}, coverage = {coverage:.2f}x\n")
                log_file.write(f"Processed {fastq}\n")

rule align_controls:
    input:
        fastq="data/{sample}.fastq",
        reference="data/reference.fasta"
    output:
        bam="results/bam/{sample}.bam"
    log:
        "results/logs/{sample}_alignment.log"
    shell:
        """
        bwa mem {input.reference} {input.fastq} |
        samtools view -Sb - > {output.bam} 2> {log}
        """

rule sort_bam_controls:
    input:
        bam="results/bam/{sample}.bam"
    output:
        sorted_bam="results/bam/{sample}_sorted.bam"
    log:
        "results/logs/{sample}_sort.log"
    shell:
        """
        samtools sort {input.bam} -o {output.sorted_bam} > {log} 2>&1
        """

rule index_sorted_bam_controls:
    input:
        sorted_bam="results/bam/{sample}_sorted.bam"
    output:
        bam_index="results/bam/{sample}_sorted.bam.bai"
    log:
        "results/logs/{sample}_index.log"
    shell:
        """
        samtools index {input.sorted_bam} > {log} 2>&1
        """

rule create_mpileup_control:
    input:
        bam="results/bam/{sample}_sorted.bam",
        bam_index="results/bam/{sample}_sorted.bam.bai",
        reference="data/reference.fasta"
    output:
        mpileup="results/mpileup/{sample}.mpileup"
    log:
        "results/logs/{sample}_mpileup.log"
    shell:
        """
        samtools mpileup -d 100000 -f {input.reference} {input.bam} > {output.mpileup} 2> {log}
        """

rule run_varscan_control:
    input:
        mpileup="results/mpileup/{sample}.mpileup"
    output:
        variants_vcf="results/vcf/{sample}_rare_variants.vcf"
    log:
        "results/logs/{sample}_varscan.log"
    shell:
        """
        java -jar /Users/anapomerash/bioinformatics_project/tools/VarScan.v2.3.9.jar mpileup2snp {input.mpileup} \
        --min-var-freq 0.001 \
        --output-vcf 1 > {output.variants_vcf} 2> {log}
        """
rule annotate_control_variants:
    input:
        vcf="results/vcf/{sample}_rare_variants.vcf"
    output:
        annotated_vcf="results/vcf/{sample}_rare_variants_annotated.vcf"
    log:
        "results/logs/{sample}_annotate.log"
    shell:
        """
        java -jar /Users/anapomerash/miniforge3/envs/myenv_x86/share/snpeff-5.2-1/snpEff.jar Flu {input.vcf} > {output.annotated_vcf} 2> {log}
        """

rule parse_vcf_to_csv:
    input:
        vcf="results/vcf/{sample}_rare_variants_annotated.vcf"
    output:
        csv="results/csv/{sample}_parsed_variants.csv"
    log:
        "results/logs/{sample}_parse_vcf.log"
    run:
        import csv
        import vcfpy

        # Открытие VCF файла
        vcf_reader = vcfpy.Reader.from_path(input.vcf)
        output_file = output.csv

        # Создаем CSV-файл с нужными колонками
        with open(output_file,'w',newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Reference", "Position", "Alternate", "Frequency", "Annotation"])

            for record in vcf_reader:
                ref = record.REF
                pos = record.POS
                alt = ",".join([alt.value for alt in record.ALT])
                freq = record.calls[0].data.get("FREQ","N/A")  # Извлекаем частоту из данных вызова
                annotations = record.INFO.get("ANN",["N/A"])  # Извлекаем аннотации

                # Записываем каждую аннотацию отдельно
                for annotation in annotations:
                    writer.writerow([ref, pos, alt, freq, annotation])

        print(f"Таблица для {input.vcf} успешно создана: {output_file}")

rule create_comparison_table:
    input:
        roommate="results/csv/roommate_parsed_variants.csv",
        controls=[
            "results/csv/SRR1705858_parsed_variants.csv",
            "results/csv/SRR1705859_parsed_variants.csv",
            "results/csv/SRR1705860_parsed_variants.csv"
        ]
    output:
        table="results/csv/comparison_table.csv",
        log="results/logs/comparison_table_stats.txt"
    log:
        "results/logs/create_comparison_table.log"
    run:
        import pandas as pd

        # Чтение данных
        roommate_data = pd.read_csv(input.roommate)
        control_data = pd.concat([pd.read_csv(file) for file in input.controls], ignore_index=True)

        # Преобразование частот в числовой формат
        roommate_data["Frequency"] = roommate_data["Frequency"].str.rstrip("%").str.replace(",", ".").astype(float) / 100
        control_data["Frequency"] = control_data["Frequency"].str.rstrip("%").str.replace(",", ".").astype(float) / 100

        # Расчёт средней и стандартного отклонения для контрольных данных (без высокочастотных вариантов)
        filtered_control_data = control_data[control_data["Frequency"] <= 0.9]
        control_mean = filtered_control_data["Frequency"].mean()
        control_std = filtered_control_data["Frequency"].std()

        # Определение порога значимости
        threshold = control_mean + 3 * control_std

        # Значимые мутации в roommate
        roommate_data["Significant"] = roommate_data["Frequency"] > threshold

        # Объединение данных roommate и control
        comparison_table = roommate_data.copy()
        comparison_table["In Control"] = comparison_table.apply(
            lambda row: control_data[
                (control_data["Position"] == row["Position"]) &
                (control_data["Reference"] == row["Reference"]) &
                (control_data["Alternate"] == row["Alternate"])
            ].any().any(),
            axis=1
        )

        # Сохранение результатов
        comparison_table.to_csv(output.table, index=False)

        # Логирование
        with open(output.log, "w") as log_file:
            log_file.write(f"Средняя частота контрольных данных: {control_mean:.6f}\n")
            log_file.write(f"Стандартное отклонение контрольных данных: {control_std:.6f}\n")
            log_file.write(f"Порог значимости: {threshold:.6f}\n")
            log_file.write(f"Количество значимых мутаций в roommate: {comparison_table['Significant'].sum()}\n")

        print("Сравнительная таблица успешно создана!")

rule analyze_epitopes:
    input:
        mutations="results/csv/significant_variants.csv"
    output:
        mutations_epitopes="results/csv/epitope_mutations.csv"
    log:
        "results/logs/analyze_epitopes.log"
    run:
        import pandas as pd

        # Загрузка данных с высокоуверенными мутациями
        mutations_data = pd.read_csv(input.mutations)

        # Определение эпитопов HA и NA
        ha_epitopes = {
            "A": [122, 124, 126] + list(range(130, 134)) + [135, 137, 138, 140] + list(range(142, 147)) + [150, 152, 168],
            "B": [128, 129] + list(range(155, 161)) + [163, 165] + list(range(186, 191)) + list(range(192, 195)) + [196, 197, 198],
            "C": list(range(44, 49)) + [50, 51, 53, 54, 273, 275, 276] + list(range(278, 281)) + [294, 297, 299, 300] + list(range(304, 313)),
            "D": [96, 102, 103, 117, 121, 167] + list(range(170, 178)) + [179, 182, 201, 203] + list(range(207, 210)) +
                 list(range(212, 220)) + list(range(226, 231)) + [238, 240, 242, 244] + list(range(246, 249)),
            "E": [57, 59, 62, 63, 67, 75, 78] + list(range(80, 84)) + list(range(86, 89)) + [91, 92, 94, 109, 260, 261, 262, 265]
        }

        na_epitopes = {
            "A": list(range(383, 388)) + list(range(389, 395)) + [396, 399, 400, 401, 403],
            "B": [197, 198, 199, 200, 221, 222],
            "C": list(range(328, 333)) + [334, 336, 338, 339] + list(range(341, 345)) + [346, 347] + list(range(357, 360)) + list(range(366, 371))
        }

        # Функция для проверки принадлежности к эпитопам
        def check_epitope(position, epitopes):
            for epitope, positions in epitopes.items():
                if position in positions:
                    return epitope
            return None

        # Анализ мутаций
        mutations_data["Epitope_HA"] = mutations_data["Position"].apply(lambda pos: check_epitope(pos, ha_epitopes))
        mutations_data["Epitope_NA"] = mutations_data["Position"].apply(lambda pos: check_epitope(pos, na_epitopes))

        # Фильтрация мутаций в эпитопах
        epitope_mutations = mutations_data[(mutations_data["Epitope_HA"].notnull()) | (mutations_data["Epitope_NA"].notnull())]

        # Сохранение результатов
        epitope_mutations.to_csv(output.mutations_epitopes, index=False)

        # Логирование
        with open(str(log), "w") as log_file:
            log_file.write(f"Обнаружено {len(epitope_mutations)} мутаций в эпитопах\n")