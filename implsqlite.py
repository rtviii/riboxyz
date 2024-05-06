
from sqlalchemy import create_engine
from sqlmodel import SQLModel, Session

from ribctl.etl.ribosome_assets import RibosomeAssets


engine = create_engine("sqlite:///ribxz.db")
SQLModel.metadata.create_all(engine)



for s in RibosomeAssets.list_all_structs()[:100]:
    print(RibosomeAssets(s))


# with Session(engine) as session:
#     session.add(hero_1)
#     session.add(hero_3)
#     session.commit()
