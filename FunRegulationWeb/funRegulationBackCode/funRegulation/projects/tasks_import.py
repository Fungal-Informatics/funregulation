import logging

from celery import shared_task
from django_celery_results.models import TaskResult

from funRegulation.task_utils import FunRegulationBaseTask

import sys
from typing import List

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GenomeApi as DatasetsGenomeApi
from ncbi.datasets.metadata.genome import get_assembly_metadata_by_bioproject_accessions

from ncbi.datasets.package import dataset

def download_files_ncbi():
    accessions: List[str] = ["GCA_003184695.1"]
    zipfile_name = "fungi_reference.zip"

    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        try:
            print("Begin download of genome data package ...")
            genome_ds_download = genome_api.download_assembly_package(
                accessions,
                include_annotation_type=["PROT_FASTA"],
                _preload_content=False,
            )

            with open(zipfile_name, "wb") as f:
                f.write(genome_ds_download.data)
            print(f"Download completed -- see {zipfile_name}")
        except DatasetsApiException as e:
            sys.exit(f"Exception when calling download_assembly_package: {e}\n")

def import_genes(import_registry):
    result = task_import_genes.apply_async((import_registry.pk,))
    import_registry.task = TaskResult.objects.get(task_id=result.task_id)
    import_registry.save()

@shared_task(bind=True, name='import_genes', base=FunRegulationBaseTask)
def task_import_genes(self, import_registry_id):
    """
    Import genes task.
    """
    logging.info('importing genes based on registry %s' % import_registry_id)