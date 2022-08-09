from multiprocessing import Process, Queue, current_process, get_context
from queue import Empty
import logging
LOGGER = logging.getLogger()
# https://pythonspeed.com/articles/python-multiprocessing
CONTEXT = get_context("spawn")

def _process_generalized(queue_in, queue_out, func, *w_args):
    LOGGER.debug("_process_generalized started")
    name = current_process().name
    LOGGER.debug("%s: starting", name)
    group = queue_in.get()
    while group != "END":
        result = func(group, *w_args)
        LOGGER.debug("%s: finished batch of %d", name, len(group))
        queue_out.put(result)
        group = queue_in.get()
    LOGGER.debug("%s: ending", name)

def parallelize(iterator, threads, func, func_out, chunksize, w_args, r_args):
    group = []
    queue_in = CONTEXT.Queue(20)
    queue_out = CONTEXT.Queue()
    procs = []
    for idx in range(threads):
        procs.append(CONTEXT.Process(
            target=_process_generalized,
            name=f"worker{idx:03}",
            args=[queue_in, queue_out, func] + w_args))
    for proc in procs:
        proc.start()
    LOGGER.debug("manager: started")
    try:
        for obj in iterator:
            if len(group) < chunksize:
                group.append(obj)
            else:
                LOGGER.debug(
                    "manager: submitting batch of %d (queue in/out: %d/%d)",
                    len(group),
                    queue_in.qsize(),
                    queue_out.qsize())
                queue_in.put(group)
                group = [obj]
                try:
                    while True:
                        result = queue_out.get(False)
                        LOGGER.debug(
                            "manager: fetched result (queue in/out: %d/%d)",
                            queue_in.qsize(),
                            queue_out.qsize())
                        func_out(result, *r_args)
                        # If we don't have many results to process try to keep
                        # them busy by submitting more
                        if queue_out.qsize() < 5:
                            break
                except Empty:
                    pass
        LOGGER.debug(
            "manager: submitting last batch of %d (queue in/out: %d/%d)",
            len(group),
            queue_in.qsize(),
            queue_out.qsize())
        queue_in.put(group)
        LOGGER.debug("last batch put")
        for _ in range(threads):
            queue_in.put("END")
        LOGGER.debug("END cmds put")
        # Gather remaining results while all workers finish up.
        #
        # Note what the multiprocessing docs say about queues:
        #
        # "...if a child process has put items on a queue (and it has not used
        # JoinableQueue.cancel_join_thread), then that process will not
        # terminate until all buffered items have been flushed to the pipe.
        #
        # This means that if you try joining that process you may get a
        # deadlock unless you are sure that all items which have been put on
        # the queue have been consumed. Similarly, if the child process is
        # non-daemonic then the parent process may hang on exit when it tries
        # to join all its non-daemonic children."
        while any(p.is_alive() for p in procs):
            try:
                result = queue_out.get(False)
                LOGGER.debug("manager: queue in/out: %d/%d", queue_in.qsize(), queue_out.qsize())
                func_out(result, *r_args)
            except Empty:
                pass
        LOGGER.debug("parallelize wrapping up")
    finally:
        for proc in procs:
            if proc.is_alive():
                LOGGER.warning("manager: worker %s still alive, terminating", proc.name)
                proc.terminate()
    LOGGER.debug("parallelize done")
