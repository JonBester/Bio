#https://blog.pythonanywhere.com/169/
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

ALLOWED_EXTENSIONS = set(['fa', 'ab1', 'fasta'])
UPLOAD_FOLDER = '/Users/jonathanbester/PythonAnywhere/Uploads'

def printer(ans):
	output = ""
	for i in range(len(ans)):
		output += ans[i] 
		output += "\n"
	return output

def fasta_leader(input_data):
	for record in SeqIO.parse(input_data, "fasta"):
		wtstring_ini = record.seq
		leader_sequence = wtstring_ini
	return input_data

def gRNA_checker(num):
	num2 = "CATS" + num
	def innerchecker(num2):
		num3 = num2[0:6]
		return num3
	num4 = innerchecker(num2)
	return num4

def ab1_decode(input_data):
    #input data must be binary but given as text
    sequence = SeqIO.read(input_data, "abi")
    sequence_ini = str(sequence.seq[0:40])
    return sequence_ini

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

def handle_fa(filename):
    fa = SeqIO.to_dict(SeqIO.parse(filename, 'fasta'))
    for record in SeqIO.parse(filename, "fasta"):
        wt_seq = str(record.seq)
    wt_truncated = wt_seq[0:50]
    return wt_truncated
    #return ', '.join([str(x) for x in fa.keys()])

def ori_checker(gRNA, fasta, ab1):
    outputs_list = {}
    for record in SeqIO.parse(fasta, "fasta"):
        wt_seq = str(record.seq)
    if gRNA in wt_seq:
    	outputs_list[0] = "Yes"
    else:
    	outputs_list[0] = "No"
    outputs_list[1] = gRNA[0:10]
    outputs_list[2] = SeqIO.read(ab1, "abi")
    return outputs_list

#def do_calculation(number1, number2):
    #return number1 + number2

#def calculate_mode(number_list):
#    try:
#        return "The mode of your list is: {}".format(statistics.mode(number_list))
#    except statistics.StatisticsError as exc:
#        return "Error calculating mode: {}".format(exc)

#the below is if you want to impliment in python terminal
#inputs = []
#while True:
#    if len(inputs) != 0:
#        print("Numbers so far:")
#        for input_value in inputs:
#            print(input_value)
#    value= input("Enter a number, or just hit return to calculate: ")
#    if value == "":
#            break
#    try:
#        inputs.append(float(value))
#    except:
#        print("{} is not a number!!")
#
#print(calculate_mode(inputs))
