from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML


# query = next(SeqIO.parse("23hed1.contigs — biopython.fa", "fasta")).seq
# result = NCBIWWW.qblast("blastn", "nt", query, hitlist_size=1)

# with open('result.xml', 'w') as res_file:
#     text = result.read()
#     res_file.write(text)

# xml = NCBIXML.parse(open("result.xml"))

# xml_list = [i.alignments for i in xml]

# for i in xml_list[0]:
#     print(i)
# print(xml_list[0][0].title)

with open(
    "23hed1_filtered_taxa.tsv", "r"
) as file:  # открываем файл, полчуенные после работы mmseq и его теста
    list_with_filtered_taxa = list()
    for i in file.readlines()[
        1:
    ]:  # убираем названия столбцов и создаем цикл для перебора файла
        if "fragment" in i or "complete far" in i or "complete close" in i:
            list_with_filtered_taxa.append(
                i.split("\t")
            )  # формируем по каждой строке список, где элемент списка - значение определенного столбца
    for i in list_with_filtered_taxa:
        i[-1] = i[-1].replace(
            "\n", ""
        )  # в последнем элементе каждого списка есть \n, убираем это
    list_with_filtered_id = list(
        [i[0] for i in list_with_filtered_taxa]
    )  # формируем список из id последовательностей, которые при тестирование не содержат "unknown"


with open("23hed1.scaffolds.fa", "r") as all_scaffolds:
    all_scaffolds_list = all_scaffolds.read().split(">")
    new_scaffolds_list = list()
    result_scaffolds_list = list()
    for i in all_scaffolds_list:
        i = i.replace("\n", "seq", 1)
        i = i.split("seq")
        i[-1] = i[-1].replace("\n", "")
        new_scaffolds_list.append(i)
    del new_scaffolds_list[0]
    for i in range(0, len(new_scaffolds_list)):
        if new_scaffolds_list[i][0] in list_with_filtered_id:
            result_scaffolds_list.append(new_scaffolds_list[i])


with open("result_filtered_scaffolds.fasta", "w") as result_result_filtered_scaffolds:
    for i in result_scaffolds_list:
        result_result_filtered_scaffolds.write(">" + i[0] + "\n" + i[1] + "\n")


with open(
    "ten_result_scaffolds_for_blast.fasta", "w"
) as result_result_filtered_scaffolds:
    for i in range(10):

        result_result_filtered_scaffolds.write(
            ">"
            + result_scaffolds_list[i][0]
            + "\n"
            + result_scaffolds_list[i][1]
            + "\n"
        )

# with open("ten_scaffold_for_blast.fasta", "w") as ten_scaffold_for_blast:
#     list_scaffold_for_blast = ''
#     for i in range(0, 10):
#         one_scaffold_for_blast = ">" + str(next(query).id) + "\n" + str(next(query).seq) + "\n"
#         list_scaffold_for_blast += one_scaffold_for_blast
#     ten_scaffold_for_blast.write(list_scaffold_for_blast)

sequence_data = open("ten_result_scaffolds_for_blast.fasta").read()


# result_blast = list()


result = NCBIWWW.qblast(
    "blastn",
    "nt",
    sequence_data,
     alignments=1,
     hitlist_size=1,
     descriptions=1,
)

# result_blast.append(result)

# print(result_blast)
with open("result.xml", "w") as result_file:
    result_file.write(result.read())
#  for i in result_blast:
#      result_text = i.read()
#      result_file.write(result_text)

with open("result_blast_output.txt", "w") as result_blast_output:
    str_result_blast_output = ""
    for count, i in enumerate(NCBIXML.parse(open("result.xml"))):
        str_result_blast_output += (
            result_scaffolds_list[count][0]
            + "\t	"
            + str(i.alignments[0]).replace("\n", ' ')
            + "\n"
        )
    result_blast_output.write(str_result_blast_output)
    print(str_result_blast_output)

# def query_0_2(query):
#     o = list()
#     p = 0
#     for i in query:
#         o.append(i.seq)
#         p += 1
#         if p == 2:
#             return o


# result_blast = list()

# def blast(query):
#     global result_blast
#     for i in query:
#         result = NCBIWWW.qblast("blastn", "nt", i, alignments=1, hitlist_size=1, format_type='text')
#         result_blast.append(result)

# blast(query_0_2(query))

# print(result_blast)

# SeqIO.write("filtered_scaffolds_fasta", filtered_scaffolds, 'fasta')
# result = NCBIWWW.qblast(
#     "blastn",
#     "nt",
#     open("result_result_filtered_scaffolds.fasta").read())
