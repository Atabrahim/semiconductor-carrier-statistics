"""Small command-line interface; emit a labeled numerical result as JSON."""

import argparse
import json
from dataclasses import asdict

from .equilibrium import solve_equilibrium


def main() -> None:
    """Parse physical inputs, solve equilibrium, and write JSON to stdout."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--temperature", type=float, default=300, help="Temperature in K")
    parser.add_argument("--donors", type=float, default=0, help="Active donors in cm^-3")
    parser.add_argument("--acceptors", type=float, default=0, help="Active acceptors in cm^-3")
    parser.add_argument(
        "--statistics", choices=["fermi-dirac", "boltzmann"], default="fermi-dirac"
    )
    parser.add_argument("--ionization", choices=["incomplete", "complete"], default="incomplete")
    args = parser.parse_args()
    try:
        result = solve_equilibrium(
            args.temperature,
            args.donors,
            args.acceptors,
            statistics=args.statistics,
            ionization=args.ionization,
        )
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps({"result_type": "NUMERICAL RESULT", **asdict(result)}, indent=2))


if __name__ == "__main__":
    main()
