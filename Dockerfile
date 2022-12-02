# base image
FROM ubuntu:22.04

# ----- SETUP -------------------
# vars inside the container
ENV DJANGO=/home/backend/api
ENV PYMOL_SOURCE=/home/backend/pymol
ENV INGRESS=/home/backend/ingress


# set work directory
RUN mkdir -p $DJANGO

# --------------------------------------
# ----- DEPENDENCIES -------------------
RUN apt-get update -y && apt-get install curl git vim \
software-properties-common apt-transport-https python3 \
build-essential python3-dev python3-pip python3-venv \
libglew-dev libpng-dev libfreetype6-dev libxml2-dev \
libmsgpack-dev python3-pyqt5.qtopengl libglm-dev libnetcdf-dev \
apt-transport-https ca-certificates libssl-dev wget  -y 



# --- NEO4J ----------------
RUN add-apt-repository -y ppa:openjdk-r/ppa
RUN apt-get update

RUN add-apt-repository -y universe

RUN wget -O neo.gpg.key https://debian.neo4j.com/neotechnology.gpg.key && apt-key add --no-tty neo.gpg.key
RUN echo 'deb https://debian.neo4j.com stable latest' | tee -a /etc/apt/sources.list.d/neo4j.list
RUN apt-get update
RUN mkdir /etc/ssl/certs/java/
RUN apt install -y --reinstall -o Dpkg::Options::="--force-confask,confnew,confmiss" --reinstall ca-certificates-java ssl-cert openssl ca-certificates
RUN apt-get install -y neo4j
# --- NEO4J ----------------


RUN pip3 install virtualenv netCDF4
COPY __pymol_source $PYMOL_SOURCE



# --------------------------------------
# ----- PYMOL_INSTALLATION -------------
# add pymol libs to path
ADD __pymol_source $PYMOL_SOURCE
ENV PYMOL_PATH="${PYMOL_SOURCE}/__pymol_lib"
ENV PYTHONPATH="${PYMOL_PATH}/modules" 
RUN mkdir -p $PYMOL_PATH
WORKDIR $PYMOL_SOURCE
RUN python3 setup.py build install --home="${PYMOL_PATH}" --install-lib="${PYMOL_PATH}/modules/" --install-scripts="${PYMOL_PATH}"


# ------------------------------------
# COPY rbxz_bend $DJANGO
COPY __ingress $INGRESS
RUN chmod +x "${INGRESS}/src/update_riboxyz.ts"
# -------------------------------------- NODE AND MODULES INSTALLATION
ENV NODE_VERSION=18.12.1
RUN curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.0/install.sh | bash
ENV NVM_DIR=/root/.nvm
RUN . "$NVM_DIR/nvm.sh" && nvm install ${NODE_VERSION}
RUN . "$NVM_DIR/nvm.sh" && nvm use v${NODE_VERSION}
RUN . "$NVM_DIR/nvm.sh" && nvm alias default v${NODE_VERSION}
ENV PATH="/root/.nvm/versions/node/v${NODE_VERSION}/bin/:${PATH}"
WORKDIR ${INGRESS}
ADD __ingress/package.json ${INGRESS}
RUN npm install --no-optional && npm cache clean --force
RUN npm install -g ts-node
# ----------------------------------------------------------------------


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

ENV NEO4J_URI="bolt://neo:7687"
ENV NEO4J_USER="neo4j"
ENV NEO4J_PASSWORD="rrr"
ENV NEO4J_CURRENTDB="neo4j"

# Be careful given that both the whole django project and the main module are called `rbxz_bend`. (hence not using the $DJANGO here)
COPY api /home/backend/api



# start server
WORKDIR ${DJANGO}
# CMD python3 manage.py runserver  0.0.0.0:8000 --noreload
CMD ["gunicorn", "--bind", ":8000", "--workers", "3", "rbxz_bend.wsgi:application","--reload"]

# ------------------------------------
# BIND THE ACTUAL STRUCTURES DATA TO SYSTEM FOLDER
# docker docs: https://docs.docker.com/storage/bind-mounts/
# The following would accept the folder defined at some host folder and mount it into ribetldata in django
# --mount type=bind,source=$HOST_FOLDER,destination=/home/backend/django/ribetldata
# verify at database startup that the seed_data is present, otherwise pull from: 


