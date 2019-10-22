import os
from flask import Flask, make_response, request, session, render_template, url_for, flash, redirect, send_from_directory
from processing import printer, ab1_decode, fasta_leader, gRNA_checker, allowed_file, handle_fa, ori_checker
from processing2 import crispr_decon
from Bio import SeqIO
from werkzeug import secure_filename


UPLOAD_FOLDER = '/Users/jonathanbester/FastaPro/Uploads'
ALLOWED_EXTENSIONS = set(['fa', 'ab1'])

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config["DEBUG"] = True
app.config["SECRET_KEY"] = "98eeaef5ee3fe79dad73c12f3fe649c5"



@app.route("/", methods=["GET", "POST"])
@app.route("/submit", methods=["GET", "POST"])
def crispr_de_page():
    if request.method == "GET":
        return '''
                <html>
                    <body>
                        <h1>CRISPR Deconvoluter: File Submission</h1>
                        <p>When a sequence read is scrabbled due to heteroallelic mutations at a CRISPR-Cas9 target site, 
                        this program can be used to deconvolute the read and determine the individual genotypes.
                        This is useful for determining if both mutations are in fact frame-shifting knockouts, or
                        just point mutations, and can inform which mutants we keep for breeding or experimentation. 
                        To run the program you need to paste your gRNA sequence (do not include PAM, 5'->3' so that 
                        the cut site is 3bp from the lefthand site). You must also provide the ab1 sequence trace file of your edited gene, and a fasta file of the wild type sequence.
                        <p>Paste your protospacer sequence 5'->3' without the PAM (ie PAM would come at the end if it were included, cutsite is 3bp from the right):
                        <form method="post" action="." enctype="multipart/form-data">
                            <p><input name="string" /></p>
                            <p>Upload a fasta file of the wt gene sequence:
                            <p><input type="file" name="fasta_file" /></p>
                            <p>Upload an ab1 or abi file of your mutated gene:
                            <p><input type="file" name="ab1_file" /></p>
                            <p>Please provide the name of your lab, institution or company:
                            <p><input name="name" /></p>
                            <p><input type="submit" name="action" value="Process your sequencing data" /></p>
                            <br>
                            <h2>Advanced parameters:</h2>
                            <p>Reverse Complement Wt Sequence?
                            <input type="checkbox" name="fasta_rv" value="true" >
                            <p>Reverse Complement Protospacer?
                            <input type="checkbox" name="guide_rv" value="true" >
                            <p>Define custom search bracket? (potential deletion size checked for; default=200)
                            <p><input name="bracket" value=200 /></p>
                            <p>Define custom freshhold? (sensativity of program to heteroallelic peaks - increase for strong signals (default=75, use ~125 for high and ~25 for low)
                            <p><input name="freshhold" value=75 /></p>
                            <p>Define custom insert length ceiling? (Maximum number of bp of random insertion the program checks for; default=3)
                            <p><input name="insert_length" value=3 /></p>
                            <p>Define custom upstream length? (Number of bp before the target site that the program searches for indels; default=30. Note: Using large values (eg 300) on low quality reads may detect false positives)
                            <p><input name="upstream_length" value=30 /></p>
                            <p>Define custom downstream length? (Number of bp after the target site that the program searches for indels; default=30. Note: Using large values (eg 300) on low quality reads may detect false positives)
                            <p><input name="downstream_length" value=30 /></p>                        
                            <p>Model length? (Number of bp the software uses for prediction, higher values give greater stringency but are slower. Provide your wt fasta template is small default should be quite sufficient for accurate assignment; default=14)
                            <p><input name="length" value=14 /></p>                        
                            </form>
                    </body>
                </html>
            '''

    if request.method == "POST":
        gRNA = request.form["string"]
        NAME = request.form["name"]
        print(str(NAME))
        filename3 = open(str(NAME), "w")
        filename3.write(str(Name))
        filename3.save(os.path.join(app.config['UPLOAD_FOLDER'], filename3))
        if request.form["action"] == "Process your sequencing data":
            upfile1 = request.files['fasta_file']
            upfile2 = request.files['ab1_file']
            if upfile1 and allowed_file(upfile1.filename):
                print(upfile1.filename)
                fasta_rv = request.form.get('fasta_rv')
                guide_rv = request.form.get('guide_rv')
                filename1 = secure_filename(upfile1.filename)
                filename2 = secure_filename(upfile2.filename)
                upfile1.save(os.path.join(app.config['UPLOAD_FOLDER'], filename1))
                upfile2.save(os.path.join(app.config['UPLOAD_FOLDER'], filename2))
                bracket = request.form["bracket"]
                freshhold = request.form["freshhold"]
                insert_length = request.form["insert_length"]
                upstream_length = request.form["upstream_length"]
                downstream_length = request.form["downstream_length"]
                length = request.form["length"]
            fasta = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
            ab1 = os.path.join(app.config['UPLOAD_FOLDER'], filename2)
            answer = crispr_decon(gRNA, fasta, ab1, fasta_rv, guide_rv, bracket, freshhold, insert_length, upstream_length, downstream_length, length)
#use multiple outputs then do most of the response formating at the html level
            if len(answer) == 15:    
                return '''
                        <html>
                            <h1>CRISPR Deconvoluter: Your Results</h1>
                            <body>
                                {A}<br>
                                {B}<br>
                                {C}<br>
                                {D}<br>
                                {E}<br>
                                {F}<br>
                                {G}<br>
                                {H}<br>
                                {I}<br>
                                {J}<br>
                                {K}<br>
                                {L}<br>
                                {M}<br>
                                {N}<br>
                                {O}<br>
                                <!-- Optional JavaScript -->
                                <!-- jQuery first, then Popper.js, then Bootstrap JS -->
                                <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
                                <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js" integrity="sha384-UO2eT0CpHqdSJQ6hJty5KVphtPhzWj9WO1clHTMGa3JDZwrnQq4sF86dIHNDz0W1" crossorigin="anonymous"></script>
                                <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js" integrity="sha384-JjSmVgyd0p3pXB1rRibZUAYoIIy6OrQ6VrjIEaFf/nJGzIxFDsf4x0xIM+B07jRM" crossorigin="anonymous"></script>
                            </body>
                        </html>
                    '''.format(A=answer[0], B=answer[1], C=answer[2], D=answer[3], E=answer[4], F=answer[5], G=answer[6], H=answer[7], I=answer[8], J=answer[9], K=answer[10], L=answer[11], M=answer[12], N=answer[13], O=answer[14])
            if len(answer) == 2:
                return '''
                        <html>
                            <h1>CRISPR Deconvoluter: Your Results</h1>
                            <body>
                                {A}<br>
                                {B}<br>
                            </body>
                        </html>
                    '''.format(A=answer[0], B=answer[1])
            if len(answer) == 1:
                return '''
                        <html>
                            <h1>CRISPR Deconvoluter: Your Results</h1>
                            <body>
                                {A}<br>
                            </body>
                        </html>
                    '''.format(A=answer[0])
        else:
            output_data = ""
            output_sequence = ""
            return '''
                    <html>
                        <body>
                            Response 1:
                            {leader} 
                            Response 2:
                            {output_sequence_sentence}
                        </body>
                    </html>
                '''.format(leader="Your gRNA sequence is: " + output_data, output_sequence_sentence = "Your number is: " + output_sequence)

@app.route("/result", methods=["GET", "POST"])
def result():
    return render_template('result.html', title='Result').format(leader="Your gRNA sequence is:" + output_data, output_sequence_sentence = "Your number is:" + output_sequence)

    
if __name__ == "__main__":
	app.run(debug=True)
