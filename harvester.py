import subprocess as sp
from pathlib import Path
import json


class Harvester:
    __slots__ = ('run_config', 'main_path', 'run_config_name'
                 'harvester_template', 'scheduler', 'scriptname')

    def __init__(self,  harvester_template: str, 
                 scheduler: str, run_config=None, scriptname='harvest.sh',
                 run_config_name='hconfig.json'):
        # Expects a dict to be fed to harvester_template.format(**run_config)
        self.run_config = run_config
        self.run_config_name = run_config_name
        # Template shell script, optionally with str.format slots for run_dict.
        self.harvester_template = harvester_template
        self.scheduler = scheduler
        self.scriptname = scriptname

    def prep_and_write_inputs(self, current_dir: Path):
        if self.run_config:
            harvest_script = self.harvester_template.format(**self.run_config)
            harvest_config_p = current_dir/self.run_config_name
            harvest_config_p.write_text(json.dumps(self.run_config))
        else:
            harvest_script = self.harvester_template
        harvest_script_p = current_dir/self.scriptname
        harvest_script_p.write_text(harvest_script)
        return harvest_script_p

    def reap(self, current_dir):
        harvest_script_p = self.prep_and_write_inputs(current_dir)
        try:
            with harvest_script_p.open() as f:
                scheduler_output = sp.check_output(
                    self.scheduler, stdin=f, cwd=current_dir, text=True
                )
        except sp.CalledProcessError as err:
            print(f'{self.scheduler} call threw error', err.stdout, err.stderr)
            raise
        return scheduler_output
