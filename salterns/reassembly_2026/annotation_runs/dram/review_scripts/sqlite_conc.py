"""Reviewer check: 8 concurrent DRAM-style readers of description_db.sqlite on NFS (read-only queries)."""
import sys, time, random, logging
from multiprocessing import Pool
from mag_annotator.database_handler import DatabaseHandler
from sqlalchemy import text

def work(seed):
    random.seed(seed)
    h = DatabaseHandler(logging.getLogger("x"))
    h.start_db_session() if not hasattr(h, "session") else None
    ids = [r[0] for r in h.session.execute(text("SELECT id FROM peptidase_description LIMIT 3000 OFFSET :o"), {"o": seed * 50000})]
    vids = [r[0] for r in h.session.execute(text("SELECT id FROM viral_description LIMIT 3000 OFFSET :o"), {"o": seed * 20000})]
    t0, n = time.time(), 0
    for _ in range(6):
        n += len(h.get_descriptions(random.sample(ids, 1500), "peptidase_description"))
        n += len(h.get_descriptions(random.sample(vids, 1500), "viral_description"))
        n += len(h.get_descriptions(["GH5", "GT2", "CBM50", "GH13"], "dbcan_description"))
    return seed, n, round(time.time() - t0, 1)

if __name__ == "__main__":
    with Pool(8) as p:
        for r in p.map(work, range(8)):
            print("reader", *r)
    print("SQLITE_CONCURRENT_OK")
