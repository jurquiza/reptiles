import time
import os
from celery import Celery
from celery.utils.log import get_task_logger

from my_functions import *

logger = get_task_logger(__name__)

app = Celery('tasks', broker='redis://redis:6379/0', backend='redis://redis:6379/0')


@app.task(bind=True)
def longtime_add(self, backbones_paths, domesticated_paths, max_length, num_sequences):

    backbones = get_backbones(backbones_paths)
    sequences = get_domesticated_fragments(domesticated_paths)

    concatenated_sequence2, combined_ids_list = find_and_concatenate_sequences(sequences, max_length, num_sequences)
    dfs = build_tree(combined_ids_list, sequences)
    csv_data = convert_to_csv(dfs)

    generate_reports(csv_data, backbones_paths, domesticated_paths)


    return len(csv_data)
