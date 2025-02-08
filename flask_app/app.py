from flask import Flask
from celery import Celery
import os
import json
import logging
import traceback
import time
import zipfile
import pandas as pd
import dnacauldron as dc
import threading
from Bio.Seq import Seq
from io import BytesIO, StringIO
from flask import Flask, jsonify, render_template, request, send_from_directory, session, make_response, redirect, url_for, send_file
from my_functions import get_backbones, find_and_concatenate_sequences, build_tree, get_domesticated_fragments, save_uploaded_file
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from flask_sqlalchemy import SQLAlchemy
from werkzeug.utils import secure_filename
import shutil

UPLOAD_FOLDER = "/shared/input"  # Use the shared volume for file storage
ZIP_FILE_PATH = "/shared/assembly_plan.zip" 
OUTPUT_FOLDER = "/shared/output" 

os.makedirs(UPLOAD_FOLDER, exist_ok=True)

# Set Flask's temporary directory to the shared volume
tempdir = UPLOAD_FOLDER  


#from joblib import Parallel, delayed
session = {}

app = Flask(__name__)
app.debug = True
celery = Celery('simple_worker', broker='redis://redis:6379/0', backend='redis://redis:6379/0')
app.secret_key = b'_5#y2L"F4Q8z\n\xec]/'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in {'gb'}

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        backbones_paths = request.files.getlist('backbones_directory')
        backbones_files = []
        domesticated_paths = request.files.getlist('domesticated_fragments_directory')
        domesticated_files = []

        for file in backbones_paths:
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(filepath)
                backbones_files.append(filepath)

        for file in domesticated_paths:
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(filepath)
                domesticated_files.append(filepath)

        session['domesticated_paths'] = domesticated_files
            

        # Convert uploaded files to a list of file paths
        #backbones_paths = [save_uploaded_file(file) for file in backbones_files]
        #domesticated_paths = [save_uploaded_file(file) for file in domesticated_files]
        
        max_length = int(request.form['max_length'])
        num_sequences = int(request.form['num_sequences'])

        sequences = get_domesticated_fragments(domesticated_files)
        backbones = get_backbones(backbones_files)
        
        concatenated_sequence2, combined_ids_list = find_and_concatenate_sequences(sequences, max_length, num_sequences)
        
        # Check if the concatenated sequence exceeds the max length
        #if len(concatenated_sequence2) > max_length:
        #    warning_message = f"Warning: The input sequences can only be assembled up to {len(concatenated_sequence2)} bp, which is less than the provided max length of {max_length} bp."
        #    return render_template('index.html', warning_message=warning_message)
        
        
        dfs = build_tree(combined_ids_list, sequences)

        dfs_json_list = [df.to_json() for df in dfs]
        session['dfs'] = dfs_json_list
        session['domesticated_paths'] = domesticated_files
        session['backbones_paths'] = backbones_files
        session['max_length'] = max_length
        session['num_sequences'] = num_sequences

        return render_template('result.html', dfs = dfs, enumerate = enumerate)

    return render_template('index.html')

@app.route('/download_csv/<int:index>')
def download_csv(index):
    """
    Download a DataFrame as a CSV file.

    Args:
        index (int): Index of the DataFrame to be downloaded from the session.

    Returns:
        flask.Response: A Flask response object with the CSV file for download.

    Raises:
        Exception: If an error occurs during the CSV file download process.
    """
    try:
        dfs_json_list = session.get('dfs')
        dfs = [pd.read_json(df_json) for df_json in dfs_json_list]

        df = dfs[index]
        # Create a BytesIO object to hold the CSV data
        csv_bytes = BytesIO()

        # Write the DataFrame to the BytesIO object as a CSV file
        df.to_csv(csv_bytes, index=False, encoding='utf-8')

        # Set the file pointer to the beginning of the BytesIO object
        csv_bytes.seek(0)

        # Create a Flask response
        response = make_response(csv_bytes.getvalue())

        # Set response headers for file download
        response.headers['Content-Type'] = 'text/csv'
        response.headers['Content-Disposition'] = f'attachment; filename=AssemblyPlan_{index}.csv'

        return response
    
    except Exception as e:
        # Log the exception for debugging
        logging.error(f"Error downloading CSV file: {e}")
        return jsonify({'error': 'Failed to download CSV file'}), 500  

@app.route('/download_all_csv')
def download_all_csv():
    """
    Download all DataFrames as a single ZIP file.

    Returns:
        flask.Response: A Flask response object with the ZIP file for download.

    Raises:
        Exception: If an error occurs during the ZIP file creation process.
    """
    try:
        dfs_json_list = session.get('dfs')
        dfs = [pd.read_json(df_json) for df_json in dfs_json_list]

        # Create a BytesIO object to hold the ZIP data
        zip_bytes = BytesIO()

        with zipfile.ZipFile(zip_bytes, 'w') as zip_file:
            for index, df in enumerate(dfs):
                # Create a CSV file in memory
                csv_buffer = StringIO()
                df.to_csv(csv_buffer, index=False)
                csv_buffer.seek(0)

                # Add the CSV file to the ZIP archive
                zip_file.writestr(f'AssemblyPlan_{index + 1}.csv', csv_buffer.getvalue())

        # Set the file pointer to the beginning of the BytesIO object
        zip_bytes.seek(0)

        # Create a Flask response
        response = make_response(zip_bytes.getvalue())

        # Set response headers for file download
        response.headers['Content-Type'] = 'application/zip'
        response.headers['Content-Disposition'] = 'attachment; filename=AssemblyPlans.zip'

        return response

    except Exception as e:
        logging.error(f"Error downloading all CSV files: {e}")
        return jsonify({'error': 'Failed to download all CSV files'}), 500



@app.route('/run-assembly-simulation', methods=['POST'])
def run_assembly():
    backbones_paths = session.get('backbones_paths', [])
    domesticated_paths = session.get('domesticated_paths', [])
    max_length = session.get('max_length', [])
    num_sequences = session.get('num_sequences', [])

    if not backbones_paths or not domesticated_paths:
        return jsonify({'error': 'No files uploaded'}), 400

    # Send file paths to Celery task
    task = celery.send_task('tasks.longtime_add', args=[backbones_paths, domesticated_paths, max_length, num_sequences])
    
    session['task_id'] = task.id  # Save task ID in session
    return jsonify({'task_id': task.id})

# Progress route
@app.route('/progress')
def progress():
    return render_template('progress.html')

@app.route('/check_progress')
def check_progress():
    task_id = session.get('task_id')
    if not task_id:
        return jsonify({'error': 'No task ID found'}), 400

    def event_stream():
        while True:
            result = celery.AsyncResult(task_id)

            if result.state == 'SUCCESS':
                yield "data: complete\n\n"
                break
            elif result.state == 'PROGRESS':
                yield f"data: {result.info.get('progress', 0)}%\n\n"
            else:
                yield "data: in_progress\n\n"

            time.sleep(1)

    return app.response_class(event_stream(), content_type='text/event-stream')


@app.route('/assembly_complete')
def assembly_complete():
    return render_template('assembly_complete.html')


@app.route('/download_assembly_plan')
def download_assembly_plan():
    """ Zips the assembly plan folder and serves it for download. """
    
    # Ensure the output folder exists
    if not os.path.exists(OUTPUT_FOLDER):
        return "Assembly plan folder not found", 404

    # Create a ZIP archive of the assembly plan
    shutil.make_archive(ZIP_FILE_PATH.replace('.zip', ''), 'zip', OUTPUT_FOLDER)

    # Serve the ZIP file for download
    return send_file(ZIP_FILE_PATH, as_attachment=True, download_name="assembly_plan.zip")


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)