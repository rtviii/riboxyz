from django.contrib.auth.models import User
from django.shortcuts import get_object_or_404
from ninja import Router
from schema import data_requests as reqs
from schema import data_responses as resps

router = Router()

@router.get('/resp/{resp_id}')
def get_article(request, resp_id: int):
    article = [ 1,2,3 ]
    return article


@router.get('/resps', response=list[RibosomeResponse])
def get_articles(request):
    # articles = Article.objects.all()
    return ["a", "b", "c"]
