"""Offline reconstruction from the immutable existing presentation source."""
from .io import *
import base64
import mimetypes
import re
from bs4 import BeautifulSoup

def build(args, out):
    root = Path(args.clinical_root)
    source = root/"presentation/source/index.html"
    html = source.read_text(encoding="utf-8")
    soup = BeautifulSoup(html, "html.parser")
    index = []
    def inline(path):
        p = (source.parent/path).resolve()
        if source.parent.resolve() not in p.parents or not p.is_file():
            raise ValueError("Missing or unsafe local presentation asset: "+path)
        mime = mimetypes.guess_type(p.name)[0] or "application/octet-stream"
        index.append(dict(path=path, bytes=p.stat().st_size, sha256=sha(p)))
        return "data:"+mime+";base64,"+base64.b64encode(p.read_bytes()).decode()
    for tag in soup.find_all(src=True):
        path = tag["src"]
        if not path.startswith("data:"):
            if "://" in path:
                raise ValueError("Offline reconstruction refuses remote presentation assets.")
            tag["src"] = inline(path)
    for style in soup.find_all("style"):
        text = style.get_text()
        text = re.sub(r'url\(["\']?([^)"\']+)["\']?\)',
            lambda m: "url("+inline(m[1])+")" if not m[1].startswith("data:") else m[0], text)
        style.string = text
    ids = [x["data-slide-id"] for x in soup.select("section[data-slide-id]")]
    if len(ids) != len(set(ids)) or len(ids) != 16:
        raise ValueError("Existing 16-slide source identity check failed.")
    for figure in root.glob("interaction_figures/*.data.tsv.gz"):
        index.append(dict(path=str(figure.relative_to(root)), bytes=figure.stat().st_size, sha256=sha(figure)))
    packed = root/"Alcohol_Microbiome_Main_Results.deck.html"
    template = BeautifulSoup(packed.read_text(encoding="utf-8"), "html.parser")
    existing = template.select("section[data-slide-id]")
    if [x["data-slide-id"] for x in existing] != ids:
        raise ValueError("Source and existing packed slide identities differ.")
    for old, new in zip(existing, soup.select("section[data-slide-id]")):
        normalize = lambda s: re.sub(r"\s+", " ", s.get_text(" ", strip=True))
        if normalize(old) != normalize(new):
            raise ValueError("Source narrative differs from the immutable packed deck.")
        images = old.find_all("img")
        replacements = new.find_all("img")
        if len(images) != len(replacements):
            raise ValueError("Slide asset structure differs from the frozen packed deck.")
        for image, replacement in zip(images, replacements):
            image["src"] = replacement["src"]
    template.find_all("style")[-1].string = soup.find("style").get_text()
    index.append(dict(path=packed.name, bytes=packed.stat().st_size, sha256=sha(packed)))
    result = str(template)
    (out/"Alcohol_Microbiome_Main_Results.deck.html").write_text(result, encoding="utf-8")
    save(pd.DataFrame(index), out/"asset_provenance.tsv")
    dump(dict(slide_ids=ids, narrative="Immutable existing source, not translated",
        standalone_PPTX="NOT_RUN", standalone_PDF="NOT_RUN",
        export_route="Open the HTML presentation in Phylo and use Export to PowerPoint or PDF."),
        out/"export_status.json")
    return dict(slides=len(ids), assets=len(index), offline=True, preserved_source_sha256=sha(source),
        preserved_editable_ids=len(template.select("[data-phylo-id]")),
        runtime="Preserved from the immutable user-provided packed presentation, not a missing external packer.")
