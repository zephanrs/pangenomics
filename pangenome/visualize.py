# visualize.py
import graphviz
import matplotlib.pyplot as plt
from .graph import Graph
from .seed import Seed
from .align import Edit

def _palette(n: int):
  cmap = plt.get_cmap("tab10")
  return [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
          for r, g, b, _ in (cmap(i) for i in range(n))]

ALIGN_COLOR = _palette(4)[3]  # tab10 red

def create_graph(graph: Graph, name: str, hedges=None):
  dot = graphviz.Digraph(filename=name, engine='dot', format='png')
  dot.attr(rankdir='LR', nodesep='0.5', ranksep='0.8')
  dot.attr('node', shape='box', fontname='Helvetica',
           style='filled', fillcolor='white')
  dot.attr(dpi='300')

  # graph nodes
  for id, seq in graph.nodes.items():
    label = seq
    dot.node(id, label=label)

  # graph edges
  for src, dsts in graph.edges.items():
    for dst in dsts:
      u, v = src.node, dst.node
      srcp = "w" if src.rev else "e"
      dstp = "e" if dst.rev else "w"

      if hedges and (u, v) in hedges:
        dot.edge(f"{u}:{srcp}", f"{v}:{dstp}",
                 color=ALIGN_COLOR, penwidth="3")
      else:
        dot.edge(f"{u}:{srcp}", f"{v}:{dstp}", color="black")

  return dot

def graph2png(graph: Graph, name: str, view: bool = False):
  dot = create_graph(graph, name)
  dot.render(cleanup=True, view=view)

def draw_seeds(dot, graph, query, seeds, k):
  # unique k-mers
  unique_kmers = []
  for s in seeds:
    kmer = query[s.qpos : s.qpos + k]
    if kmer not in unique_kmers:
      unique_kmers.append(kmer)
  
  # limit to 8
  palette = _palette(min(8, len(unique_kmers)))
  color_map = {kmer: palette[i] for i, kmer in enumerate(unique_kmers[:8])}

  # group highlighting
  node_spans = {}
  for seed in seeds:
    kmer = query[seed.qpos : seed.qpos + k]
    if kmer in color_map:
      color = color_map[kmer]
      node_spans.setdefault(seed.node, []).append((seed.npos, seed.npos + k, color))

  for node, spans in node_spans.items():
    seq = graph.nodes[node]
    chars = list(seq)
    for a, b, color in spans:
      for i in range(a, b):
        chars[i] = f'<U><FONT COLOR="{color}">{chars[i]}</FONT></U>'
    dot.node(node, label=f'<<FONT FACE="Helvetica">{"".join(chars)}</FONT>>')

  # underline
  q = list(query)
  for seed in seeds:
    kmer = query[seed.qpos : seed.qpos + k]
    if kmer in color_map:
      color = color_map[kmer]
      for i in range(seed.qpos, seed.qpos + k):
        q[i] = f'<U><FONT COLOR="{color}">{q[i]}</FONT></U>'

  # query node
  dot.node(
    "query",
    label=(
      '<<TABLE ALIGN="CENTER" BORDER="0" CELLBORDER="0">'
      f'<TR><TD><FONT FACE="Helvetica" POINT-SIZE="20">{"".join(q)}</FONT></TD></TR>'
      '</TABLE>>'
    ),
    shape="plaintext"
  )

  with dot.subgraph() as s:
    s.attr(rank="sink")
    s.node("query")
    
def seeds2png(graph: Graph, query: str, seeds: list[Seed],
              k: int, name: str, view=False):
  dot = create_graph(graph, name)
  draw_seeds(dot, graph, query, seeds, k)
  dot.render(cleanup=True, view=view)

def align2png(graph: Graph, query: str, edits: list[Edit],
              name: str, view=False):

  # ordered path of nodes
  path = []
  for e in edits:
    if not path or path[-1] != e.node:
      path.append(e.node)

  hedges = set(zip(path, path[1:]))

  dot = create_graph(graph, name, hedges=hedges)

  # recolor existing nodes
  for node in path:
    dot.node(node, color=ALIGN_COLOR, penwidth="3")

  dot.render(cleanup=True, view=view)

def alignment_strings(graph: Graph, query: str, edits: list[Edit]):
  q = []
  r = []

  for e in edits:
    if e.op in ("M", "X"):
      q.append(query[e.qpos])
      r.append(graph.nodes[e.node][e.npos])
    elif e.op == "I":
      q.append(query[e.qpos])
      r.append("-")
    elif e.op == "D":
      q.append("-")
      r.append(graph.nodes[e.node][e.npos])

  return "".join(q), "".join(r)
