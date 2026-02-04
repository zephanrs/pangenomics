import cmd, os, glob
from .graph import gfa2graph, graph2gfa, Graph
from .seed import Seeder
from .filter import Filterer
from .align import Aligner
from .pangenome import add_alignment
from .visualize import seeds2png, align2png, graph2png, alignment_strings

class PangenomeShell(cmd.Cmd):
    intro = """
    Pangenome CLI. Available commands:
    ----------------------------------
    load <path>   : Load GFA file
    store <path>  : Save GFA file
    graph         : View current graph
    seed <seq>    : Set query and find seeds (view img)
    filter        : Filter current seeds (view img)
    extend        : Align using filtered seeds (view img + strings)
    add <name>    : Add current alignment to graph (view img)
    params        : Show parameters
    k <val>       : Set k-mer size
    w <val>       : Set window size
    reward <val>  : Set alignment reward
    help          : Show this help message
    quit          : Exit and cleanup
    """
    prompt = "(pangenome) "
    
    def __init__(self):
        super().__init__()
        self.g = Graph()
        self.query = ""
        self.seeds = []
        self.edits = []
        self.k = 3
        self.w = 5
        self.reward = 10
        self.tmp = "temp_view"

    def cleanup(self):
        for f in glob.glob(f"{self.tmp}*"):
            try: os.remove(f)
            except: pass
        print("Cleaned up temporary files.")

    def do_help(self, arg):
        print(self.intro)

    def do_load(self, arg):
        try:
            self.g = gfa2graph(arg)
            print(f"Loaded {len(self.g.nodes)} nodes.")
        except Exception as e: print(f"Error: {e}")

    def do_store(self, arg):
        try:
            graph2gfa(self.g, arg)
            print(f"Saved to {arg}")
        except Exception as e: print(f"Error: {e}")

    def do_graph(self, arg):
        graph2png(self.g, self.tmp, view=True)

    def do_seed(self, arg):
        if not arg: return print("Usage: seed <sequence>")
        self.query = arg
        seeder = Seeder(self.g, self.k, self.w)
        self.seeds = seeder.seed(self.query)
        print(f"Found {len(self.seeds)} seeds.")
        seeds2png(self.g, self.query, self.seeds, self.k, self.tmp, view=True)

    def do_filter(self, arg):
        if not self.seeds: return print("Run 'seed' first.")
        filt = Filterer(self.g, "max", self.reward, {})
        self.seeds = filt.filter(self.seeds)
        print(f"Filtered to {len(self.seeds)} seeds.")
        seeds2png(self.g, self.query, self.seeds, self.k, self.tmp, view=True)

    def do_extend(self, arg):
        if not self.seeds: return print("Run 'filter' first.")
        aligner = Aligner(self.g)
        self.edits = aligner.align(self.query, self.seeds)
        print(f"Aligned: {len(self.edits)} edits.")
        
        q_str, r_str = alignment_strings(self.g, self.query, self.edits)
        print(f"Q: {q_str}\nR: {r_str}")
        align2png(self.g, self.query, self.edits, self.tmp, view=True)

    def do_add(self, arg):
        if not self.edits: return print("Run 'extend' first.")
        name = arg if arg else "read"
        self.g = add_alignment(self.g, self.query, self.edits, name)
        print(f"Alignment '{name}' added. Graph nodes: {len(self.g.nodes)}")
        graph2png(self.g, self.tmp, view=True)

    def do_params(self, arg):
        print(f"k={self.k}, w={self.w}, reward={self.reward}")

    def do_k(self, arg): self.k = int(arg); print(f"k set to {self.k}")
    def do_w(self, arg): self.w = int(arg); print(f"w set to {self.w}")
    def do_reward(self, arg): self.reward = int(arg); print(f"reward set to {self.reward}")
    
    def do_quit(self, arg):
        self.cleanup()
        return True

if __name__ == '__main__':
    shell = PangenomeShell()
    try:
        shell.cmdloop()
    except KeyboardInterrupt:
        print("\nInterrupted.")
        shell.cleanup()