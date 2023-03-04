from django.contrib.auth.models import User
from django.shortcuts import get_object_or_404

from ninja import Router


router = Router()

@router.get('/resp/{resp_id}')
def get_article(request, resp_id: int):
    article = [ 1,2,3 ]
    return article


@router.get('/resps', response=list)
def get_articles(request):
    return ["a", "b", "c"]
