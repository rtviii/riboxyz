# base image
FROM ubuntu:22.04

# ----- SETUP -------------------
# vars inside the container
ENV DJANGO=/home/backend/api
ENV PYMOL_SOURCE=/home/backend/pymol
ENV CLI=/home/backend/ingress

# set work directory
RUN mkdir -p $DJANGO

# ----- LIBS -------------------------------------------
RUN apt-get update -y && apt-get install curl git vim \
software-properties-common apt-transport-https python3 \
build-essential python3-dev python3-pip python3-venv \
libglew-dev libpng-dev libfreetype6-dev libxml2-dev \
libmsgpack-dev python3-pyqt5.qtopengl libglm-dev libnetcdf-dev \
apt-transport-https ca-certificates libssl-dev wget  -y 
# -------------------------------------------------------------


# --- NEO4J --------------------------------------------------------------------------------
RUN add-apt-repository -y ppa:openjdk-r/ppa
RUN apt-get update

RUN add-apt-repository -y universe

RUN wget -O neo.gpg.key https://debian.neo4j.com/neotechnology.gpg.key && apt-key add --no-tty neo.gpg.key
RUN echo 'deb https://debian.neo4j.com stable latest' | tee -a /etc/apt/sources.list.d/neo4j.list
RUN apt-get update
RUN mkdir /etc/ssl/certs/java/
RUN apt install -y --reinstall -o Dpkg::Options::="--force-confask,confnew,confmiss" --reinstall ca-certificates-java ssl-cert openssl ca-certificates
RUN apt-get install -y neo4j
# --- NEO4J ---------------------------------------------------------------------------------

# ----- PYMOL_INSTALLATION -------------
RUN pip3 install virtualenv netCDF4
COPY pymol_source $PYMOL_SOURCE
ADD pymol_source $PYMOL_SOURCE
ENV PYMOL_PATH="${PYMOL_SOURCE}/__pymol_lib"
ENV PYTHONPATH="${PYMOL_PATH}/modules" 
RUN mkdir -p $PYMOL_PATH
WORKDIR $PYMOL_SOURCE
RUN python3 setup.py build install --home="${PYMOL_PATH}" --install-lib="${PYMOL_PATH}/modules/" --install-scripts="${PYMOL_PATH}"
# --- NEO4J ---------------------------------------------------------------------------------




# WORKDIR $DJANGO
ADD api/reqs.txt $DJANGO
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip3 install -r ${DJANGO}/reqs.txt
RUN pip3 install pytz   --upgrade
RUN pip3 install tzdata --upgrade
RUN pip3 install gunicorn

# port where the Django app runs
EXPOSE 8000 8001 8002

ENV NEO4J_URI       = "bolt://neo:7687"
ENV NEO4J_USER      = "neo4j"
ENV NEO4J_PASSWORD  = "rrr"
ENV NEO4J_CURRENTDB = "neo4j"

# -------------------------------------- NODE AND MODULES INSTALLATION
COPY cli $CLI
# RUN chmod +x "${INGRESS}/src/update_riboxyz.ts"

ENV NODE_VERSION=18.12.1
RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.0/install.sh | bash
ENV NVM_DIR=/root/.nvm
RUN . "$NVM_DIR/nvm.sh" && nvm install ${NODE_VERSION}
RUN . "$NVM_DIR/nvm.sh" && nvm use v${NODE_VERSION}
RUN . "$NVM_DIR/nvm.sh" && nvm alias default v${NODE_VERSION}
ENV PATH="/root/.nvm/versions/node/v${NODE_VERSION}/bin/:${PATH}"
WORKDIR ${CLI}
ADD cli/package.json ${CLI}
RUN npm install --no-optional && npm cache clean --force
RUN npm install -g ts-node
# ----------------------------------------------------------------------

# Be careful given that both the whole django project and the main module are called `rbxz_bend`. (hence not using the $DJANGO here)
COPY api /home/backend/api

# start server
WORKDIR ${DJANGO}
CMD ["gunicorn", "--bind", ":8000", "--workers", "3", "rbxz_bend.wsgi:application","--reload"]