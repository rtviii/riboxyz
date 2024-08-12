
# DevOps:

- one-click docker deployment
- how can we benefit from github actions?
- automatic update per struct -- cron job every 24 weeks



### Pymol

Pymol, so long as it serves its purpose for filling in the gaps of visual and geometric processing, should be installed the following way (taken verbatim from `docker-compose` file):

```bash
COPY pymol_source $PYMOL_SOURCE
ADD pymol_source $PYMOL_SOURCE
ENV PYMOL_PATH="${PYMOL_SOURCE}/__pymol_lib"
ENV PYTHONPATH="${PYMOL_PATH}/modules" 
RUN mkdir -p $PYMOL_PATH

WORKDIR $PYMOL_SOURCE
RUN python3 setup.py build install --home="${PYMOL_PATH}" --install-lib="${PYMOL_PATH}/modules/" --install-scripts="${PYMOL_PATH}"
# somewhere in the ~/.XXXrc's:
export PYTHONPATH="${PYTHONPATH}:${PYMOL_PATH}/modules"
```
Looking at the above by parts:
- OSS pymol installation requires the mmtf cpp packages. We grab those as is and they are already inside `pymol_source/include`.

- we install pymol headlessly throguh `setup.py` specifying a particular $PYMOL_PATH (a path of my choosing, it's only important that this ends up on pythonpath)

- to use pymol as a library, what we are really interested in is `$PYMOL_PATH/module` dir. We add that to $PYTHONPATH and voila, `from pymol import cmd` works anywhere.


The following pathing is elaborated on here: https://sourceforge.net/p/pymol/mailman/message/35988916/