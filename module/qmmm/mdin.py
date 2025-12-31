from pathlib import Path
from typing import List, Dict
from functools import lru_cache

class MdinGenerator:
    __slots__ = ('_templates',)

    def __init__(self):
        self._templates = {}
        self._init_templates()

    def _init_templates(self):
        self._templates['min'] = """Minmize all the hydrogens
 &cntrl
 ifqnt=1,
 imin=1,
 maxcyc=1,
 ncyc=1,
 ntb=1,
 ntp=0,
 ntf=1,
 ntwx= 1,
 ntpr= 1,
 ntwr= 1,
 cut=  10.0,
 ntr=1,
 restraintmask=":1-9999 & !@H=",
 restraint_wt=1.0,
 ntxo=1,
 ioutfm=0,
 /
 &qmmm
 qmmask = ':111,55',
 qmcharge = 0,
 qm_theory = 'XTB',
 spin = 2,
 verbosity = 2,
 printcharges = 1,
 qmcut = 10.0,
 qm_ewald = 1,
 writepdb = 1,
 /
 &xtb
 qm_level = 'GFN2-xTB',
 tfermi = 300.0,
 accuracy = 1.d-3
 maxiter = 500
 mmhardness = 0.0
 debug = F
 /
"""
        
        self._templates['sp'] = """SP calculation
 &cntrl
 ifqnt=1,
 imin=6,
 nstlim=1,
 dt=0.0,
 ntb=0,
 ntp=0,
 ntf=1,
 ntwx= 1,
 ntpr= 1,
 ntwr= 1,
 cut=  10.0,
 ntr=1,
 restraintmask=":1-9999 & !@H=",
 restraint_wt=1.0,
 ntxo=1,
 ioutfm=0,
 /
 &qmmm
 qmmask = ':111,55',
 qmcharge = ,
 qm_theory = 'XTB',
 spin = 2,
 verbosity = 2,
 printcharges = 1,
 qmcut = 10.0,
 qm_ewald = 1,
 writepdb = 1,
 /
 &xtb
 qm_level = 'GFN2-xTB',
 tfermi = 300.0,
 accuracy = 1.d-3
 maxiter = 500
 mmhardness = 0.0
 debug = F
 /
"""
        
        self._templates['vdwm'] = """MMSP calculation
 &cntrl
 imin=6,
 nstlim=1,
 dt=0.0,
 ntb=0,
 ntp=0,
 ntf=1,
 ntwx= 1,
 ntpr= 1,
 ntwr= 1,
 cut=  10.0,
 ntr=1,
 restraintmask=":1-9999 & !@H=",
 restraint_wt=1.0,
 ntxo=1,
 ioutfm=0,
 /
"""
        
        self._templates['min_sys2'] = """Min SYS2 calculation
 &cntrl
 ifqnt=1,
 imin=1,
 maxcyc=1,
 ncyc=1,
 ntb=1,
 ntp=0,
 ntf=1,
 ntwx= 1,
 ntpr= 1,
 ntwr= 1,
 cut=  10.0,
 ntr=1,
 restraintmask=":1-9999 & !@H=",
 restraint_wt=1.0,
 ntxo=1,
 ioutfm=0,
 /
 &qmmm
 qmmask = ':55',
 qmcharge = 0,
 qm_theory = 'XTB',
 spin = 2,
 verbosity = 2,
 printcharges = 1,
 qmcut = 10.0,
 qm_ewald = 1,
 /
 &xtb
 qm_level = 'GFN2-xTB',
 tfermi = 300.0,
 accuracy = 1.d-3
 maxiter = 500
 mmhardness = 0.0
 debug = F
 /
"""
        
        self._templates['sp_sys2'] = """SP SYS2 calculation
 &cntrl
 ifqnt=1,
 imin=6,
 nstlim=1,
 dt=0.0,
 ntb=0,
 ntp=0,
 ntf=1,
 ntwx= 1,
 ntpr= 1,
 ntwr= 1,
 cut=  10.0,
 ntr=1,
 restraintmask=":1-9999 & !@H=",
 restraint_wt=1.0,
 ntxo=1,
 ioutfm=0,
 /
 &qmmm
 qmmask = ':55',
 qmcharge = ,
 qm_theory = 'XTB',
 spin = 2,
 verbosity = 2,
 printcharges = 1,
 qmcut = 10.0,
 qm_ewald = 1,
 /
 &xtb
 qm_level = 'GFN2-xTB',
 tfermi = 300.0,
 accuracy = 1.d-3
 maxiter = 500
 mmhardness = 0.0
 debug = F
 /
"""

    @lru_cache(maxsize=8)
    def get_template(self, template_type: str) -> str:
        return self._templates.get(template_type, self._templates['min'])

    def generate_mdin(self, output_path: Path, template_type: str='min') -> str:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(self.get_template(template_type))
        return str(output_path)

    def _generate_files_for_dir(self, output_dir: Path) -> Dict:
        output_dir = Path(output_dir)
        file_specs = [('min.mdin', 'min'), ('min_sys2.mdin', 'min_sys2'), ('sp.mdin', 'sp'), ('vdwm.mdin', 'vdwm'), ('sp_sys2.mdin', 'sp_sys2')]
        try:
            files = {f'{typ}_mdin_file': self.generate_mdin(output_dir / fname, typ) for fname, typ in file_specs}
            return {'output_dir': str(output_dir), **files, 'status': 'success'}
        except Exception as e:
            return {'output_dir': str(output_dir), 'status': 'error', 'error_message': str(e)}

    def generate_batch(self, output_dirs: List[Path]) -> List[Dict]:
        return list(map(self._generate_files_for_dir, output_dirs))

def main():
    generator = MdinGenerator()
if __name__ == '__main__':
    main()