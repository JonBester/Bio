import os
from flask import Flask, make_response, request, session, render_template, url_for, flash, redirect, send_from_directory
from processing import ab1_decode, fasta_leader, gRNA_checker, allowed_file, handle_fa 
from Bio import SeqIO
from werkzeug import secure_filename


UPLOAD_FOLDER = '/Users/jonathanbester/FastaPro/Uploads'
ALLOWED_EXTENSIONS = set(['fa', 'ab1'])

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@app.route('/', methods=['GET', 'POST'])
def upload_file():
    print('ok')
    print(dir(request))
    print("Files:", request.files)

    if request.method == 'POST':
        upfile = request.files['file']
        if upfile and allowed_file(upfile.filename):
            print(upfile.filename)
            filename = secure_filename(upfile.filename)
            upfile.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            #sequence = SeqIO.read(open(os.path.join(app.config['UPLOAD_FOLDER'], filename)), "fasta")
            #sequence_ini = str(sequence.seq[0:40])
            #return sequence_ini            
            return handle_fa(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            #redirect(url_for('upload_file'))

    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <p>Paste your protospacer sequence 5'->3' without the PAM (ie PA
    <form action="" method=post enctype=multipart/form-data>
      <p><input type=file name=file>
         <input type=submit value=Upload>
    </form>
    '''



if __name__ == "__main__":
   app.run(port=5000, debug=True)