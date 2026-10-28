"""Exercise native fragment requests while retaining AppTest's session/tree.

AppTest creates a new local runner for every run. Retain the fragment storage
and forward-message queue just as AppSession does, then submit a real
fragment-targeted RerunData request. Application renderers are never bypassed.
"""
from dataclasses import replace
from unittest.mock import patch

from streamlit.runtime.scriptrunner_utils.script_requests import ScriptRequests
from streamlit.testing.v1.local_script_runner import LocalScriptRunner


class FragmentSession:
    def __enter__(self):
        self.runners = []
        self.fragment_id = None
        self.patch = patch('streamlit.testing.v1.app_test.LocalScriptRunner', side_effect=self._runner)
        self.patch.start()
        return self

    def __exit__(self, *exc):
        self.patch.stop()

    def _runner(self, *args, **kwargs):
        runner = LocalScriptRunner(*args, **kwargs)
        if self.fragment_id is not None:
            previous = self.runners[-1]
            runner._fragment_storage = previous._fragment_storage
            runner._requests = ScriptRequests()
            for message in previous.forward_msgs():
                runner.forward_msg_queue.enqueue(message)
            request = runner.request_rerun
            fragment_id = self.fragment_id
            runner.request_rerun = lambda data: request(replace(data, fragment_id=fragment_id))
        self.runners.append(runner)
        return runner

    def run(self, app, widget_key):
        state = app.session_state._state
        metadata = state._get_widget_metadata(state._get_widget_id(widget_key))
        assert metadata.fragment_id, f'{widget_key} is not owned by a fragment'
        self.fragment_id = metadata.fragment_id
        try:
            return app.run()
        finally:
            self.fragment_id = None
