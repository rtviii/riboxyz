from pydantic import BaseModel


class Organism(BaseModel): 
      domain             : tuple[int, str]
      species            : tuple[int, str]
      strain             : tuple[int, str]




