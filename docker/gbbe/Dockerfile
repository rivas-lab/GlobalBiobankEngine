# Do not use ‘latest’ tag, as ‘latest’ is a moving target.
FROM python:2.7.18
# Copy in current dir  and place it under /biobankengine
COPY . /biobankengine
RUN apt-get update && apt-get upgrade -y \
	&& apt-get install -y --no-install-recommends \
	libreadline-dev \
	build-essential \
	r-base \
	&& rm -rf /var/lib/apt/lists/*
# Install the SciDB prereqs
RUN pip install -r /biobankengine/requirements.txt \
	&& cd /biobankengine/pkg_py \
	&& pip install -r requirements.txt \
	&& pip install .

COPY entrypoint.sh /entrypoint.sh
ENTRYPOINT /entrypoint.sh
