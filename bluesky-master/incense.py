import incense
from incense import ExperimentLoader

loader = ExperimentLoader(
    mongo_uri="localhost:27017", 
    db_name='db_prueba'
)

loader