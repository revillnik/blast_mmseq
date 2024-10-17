from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# name = input('Впишите название-идентификатор для файлов, полученных в ходе скрипта (оно может быть любое)')
# filtered_taxa = input(
#     "введите путь к файлу с результирующей классификацией (полученной после обработки mmseq), пишите также тип файла через '.' "
# )

# scaffolds = input("введите путь к файлу, содержащему все скаффолды, пишите также тип файла через '.'")

name = "test"
filtered_taxa = "23hed1_filtered_taxa.tsv"
scaffolds = "23hed1.scaffolds.fa"
"bat21_132_filtered_taxa.tsv"
"bat21_132.scaffolds.fa"
"23hed1.scaffolds.fa"
"23hed1_filtered_taxa.tsv"

with open(
    filtered_taxa, "r"
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
    )  # формируем список из id последовательностей, которые при тестировании показали значения после тестирования "fragment" или "complete far" или "complete close"


with open(scaffolds, "r") as all_scaffolds:  # открываем файл со всеми скаффолдами
    all_scaffolds_list = all_scaffolds.read().split(
        ">"
    )  # получаем все скаффолды в виде списка (элементы содержат строку с id + seq)
    new_scaffolds_list = list()
    result_scaffolds_list = list()
    for i in all_scaffolds_list:
        i = i.replace(
            "\n", "seq", 1
        )  # убираем первый \n, меняем на seq, чтобы мы могли разделить id и seq
        i = i.split("seq")  # делим каждый элемент списка, на 2 - id и seq
        i[-1] = i[-1].replace("\n", "")  # убираем лишние \n из seq
        new_scaffolds_list.append(i)
    del new_scaffolds_list[0]  # удаляем пустую строку (не знаю почему она тут)
    for i in range(0, len(new_scaffolds_list)):
        if (
            new_scaffolds_list[i][0] in list_with_filtered_id
        ):  # получаем список скаффолдов, но только тех, что есть в отфильтрованном списке таксонов (то есть убираем лишние скаффолды)
            result_scaffolds_list.append(new_scaffolds_list[i])


with open(
    f"{name}_result_filtered_scaffolds.fasta", "w"
) as result_result_filtered_scaffolds:  # создаем файл с расширением fasta, где будут лежать все отфильтрованные скаффолды
    for i in result_scaffolds_list:
        result_result_filtered_scaffolds.write(">" + i[0] + "\n" + i[1] + "\n")


with open(
    f"{name}_ten_result_scaffolds_for_blast.fasta", "w"
) as result_result_filtered_scaffolds:  # создаем файл с 10-ю скаффолдами из отфильтрованных (так как больше команда qblast не дает вписать)
    str_for_filtered_scaffolds = ""
    for i in range(0, len(result_scaffolds_list)):
        if len(str_for_filtered_scaffolds) > 195000:
            break
        else:
            str_for_filtered_scaffolds += (
                ">"
                + result_scaffolds_list[i][0]
                + "\n"
                + result_scaffolds_list[i][1]
                + "\n"
            )
    result_result_filtered_scaffolds.write(str_for_filtered_scaffolds)

sequence_data = open(
    f"{name}_ten_result_scaffolds_for_blast.fasta"
).read()  # получаем сырые данные по 10-ти отфильтрованным скаффолдам для blast

result = NCBIWWW.qblast(
    "blastn",
    "nt",
    sequence_data,
    megablast=True,
) # проводим blast по скаффолдам из отфильтрованного списка, при этом ограничимся 1 выравниванием (почему-то при мегабласте проходит сразу много выравниваний)

with open(
    f"{name}_result.xml", "w"
) as result_file:  # записываем данные полученные после blast в файл
    result_file.write(result.read())

with open(
    f"{name}_result_blast_output.txt", "w"
) as result_blast_output:  # создаем текстовый файл, где будет лежать краткая информация по проведенным выравниваниям
    str_result_blast_output = ""
    for count, i in enumerate(NCBIXML.parse(open(f"{name}_result.xml"))):
        try:
            str_result_blast_output += (
         result_scaffolds_list[count][0]
         + "\t	"
         + str(i.alignments[0]).replace("\n", " ")
         + "\n"
     ) # заметь, что тут мы берем первое выравнивание! (т.к. у него больше всего оценка)
        except:
           pass
    result_blast_output.write(str_result_blast_output)
