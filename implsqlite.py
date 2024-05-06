import sqlite3 
connection = sqlite3.connect("ribxz.db")
print(connection.total_changes)
cursor = connection.cursor()
cursor.execute("CREATE TABLE IF NOT EXISTS example (id INTEGER, name TEXT, age INTEGER)")
