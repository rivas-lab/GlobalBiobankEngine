import boto3

import numpy
import sklearn
import scipy

s3_client = boto3.client('s3')

def process(event, context):
    """
    Process a file upload.
    """

    # Get the uploaded file's information
    bucket = event['Records'][0]['s3']['bucket']['name'] # Will be `my-bucket`
    key = event['Records'][0]['s3']['object']['key'] # Will be the file path of whatever file was uploaded.

    # Get the bytes from S3
    #s3_client.download_file(bucket, key, '/tmp/' + key) # Download this file to writable tmp space.
    #file_bytes = open('/tmp/' + key).read()
    print numpy.version.version
    print scipy.version.version
    print sklearn.__version__
