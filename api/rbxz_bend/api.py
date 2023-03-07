from django.contrib.auth.models import User
from django.shortcuts import get_object_or_404
from ninja import Router
from ninja import Query, Schema
from ribctl.db.data import QueryOps
from schema.v0 import NeoStruct

router = Router()



@router.get('/v0/get_all_structures', response=list[NeoStruct])
def get_all_structures(request):
    qo = QueryOps()
    qo.get_all_structures()

    


    




# @router.get('/resp/{resp_id}')
# def get_article(request, resp_id: int):
#     article = [ 1,2,3 ]
#     return article


# @router.get('/resps', response=list[RibosomeResponse])
# def get_articles(request):
#     # articles = Article.objects.all()
#     return ["a", "b", "c"]
