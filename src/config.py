import os
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class RuntimeConfig:
    query_virus_dir: str
    output_dir: str
    intermediate_dir: str
    data_dir: str
    genome_list: Optional[str]
    short_contig: bool
    topN: int
    num_threads: int

    @classmethod
    def from_args(cls, args):
        query_virus_dir = os.path.abspath(os.path.expanduser(args.query_virus_dir))
        output_dir = os.path.abspath(os.path.expanduser(args.output_dir))
        intermediate_dir = os.path.abspath(os.path.expanduser(args.intermediate_dir))
        data_dir = os.path.abspath(os.path.expanduser(args.data_dir))
        genome_list = args.genome_list
        if genome_list is not None:
            genome_list = os.path.abspath(os.path.expanduser(genome_list))
        return cls(
            query_virus_dir=query_virus_dir,
            output_dir=output_dir,
            intermediate_dir=intermediate_dir,
            data_dir=data_dir,
            genome_list=genome_list,
            short_contig=args.short_contig,
            topN=args.topN,
            num_threads=args.num_Threads,
        )
