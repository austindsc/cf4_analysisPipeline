import sqlite3 as sql
import pandas as pd

import subprocess
bashCommand = "rm test_runLists.db"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

if __name__ == '__main__':
  conn = sql.connect("test_runLists.db")
  c = conn.cursor()

  c.execute("""
                  CREATE TABLE IF NOT EXISTS allRuns
                  ([runID] INTEGER PRIMARY KEY, [gammaSource] STRING, [gammaPos] FLOAT, [alpha] FLOAT, [bkg] BOOLEAN, [])
                  """)
  conn.commit()
  # runID is date and time
  # gamma source is string
  # gamms pos is float?

  # alpha is float. if non zeror its position

  c.execute('''
          INSERT INTO allRuns (runID, gamma, bkg)
                VALUES
                (202204141147, 1, 1)
          ''')
  conn.commit()

  c.execute('''
          SELECT *
          FROM allRuns;
          ''')

  print(c.fetchall())
  # df = pd.DataFrame(c.fetchall(), columns=['runID','gamma'])
  # print (df)  


