import subprocess
import os

def calculate_rnadistance(s1, s2, method='default'):
    """
    Calculates the distance between two RNA structures using the external
    ViennaRNA RNAdistance executable.

    Parameters:
    -----------
    s1, s2 : str
        The dot-bracket structure strings.
    method : str
        'default' - Tree Edit Distance (Default)
        'string'  - String Alignment (-D F)
        'bp'      - Base Pair Distance (-D p)

    Returns:
    --------
    float/int : The calculated distance.
    """
    
    exe_path = r"C:\Program Files (x86)\ViennaRNA Package\RNAdistance.exe"
    
    if not os.path.exists(exe_path):
        print(f"Error: RNAdistance not found at {exe_path}")
        return None

    cmd_args = [exe_path]
    
    if method == 'string':
        cmd_args.extend(["-D", "F"])
    elif method == 'bp':
        cmd_args.extend(["-D", "p"])
    elif method == 'default':
        pass 

    input_str = f"{s1}\n{s2}\n"

    try:
        res = subprocess.run(
            cmd_args,
            input=input_str,
            text=True,
            capture_output=True
        )
        
        output = res.stdout.strip()
        
        if not output:
            return None

        last_line = output.splitlines()[-1]
        
        if ":" in last_line:
            distance = float(last_line.split(":")[-1].strip())
        else:
            parts = last_line.split()
            if parts and parts[-1].replace('.', '', 1).isdigit():
                distance = float(parts[-1])
            else:
                print(f"Warning: Could not parse output: {last_line}")
                return None

        return int(distance)

    except Exception as e:
        print(f"System Error running RNAdistance: {e}")
        return None