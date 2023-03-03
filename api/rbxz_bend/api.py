from django.contrib.auth.models import User
from django.shortcuts import get_object_or_404
from ninja import Router
from ribctl.lib.types.types_ribosome import RibosomeResponse

router = Router()

# @router.post('/utils/resp')
# def create_article(request, payload: ArticleIn):
#     data = payload.dict()
#     try:
#         author = User.objects.get(id=data['author'])
#         del data['author']
#         article = Article.objects.create(author=author, **data)
#         return {
#             'detail': 'Article has been successfully created.',
#             'id': article.id,
#         }
#     except User.DoesNotExist:
#         return {'detail': 'The specific user cannot be found.'}


@router.get('/resp/{resp_id}')
def get_article(request, resp_id: int):
    article = [ 1,2,3 ]
    return article


@router.get('/resps', response=list[RibosomeResponse])
def get_articles(request):
    # articles = Article.objects.all()
    return ["a", "b", "c"]
