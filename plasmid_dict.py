
# ------------------ Parsing Plasmid dictionary --------------------
from typing import Optional, Any, List, TypeVar, Callable, Type, cast
from enum import Enum


T = TypeVar("T")
EnumT = TypeVar("EnumT", bound=Enum)


def from_str(x: Any) -> str:
    assert isinstance(x, str)
    return x


def from_none(x: Any) -> Any:
    assert x is None
    return x


def from_union(fs, x):
    for f in fs:
        try:
            return f(x)
        except:
            pass
    assert False


def from_int(x: Any) -> int:
    assert isinstance(x, int) and not isinstance(x, bool)
    return x


def from_list(f: Callable[[Any], T], x: Any) -> List[T]:
    assert isinstance(x, list)
    return [f(y) for y in x]


def to_enum(c: Type[EnumT], x: Any) -> EnumT:
    assert isinstance(x, c)
    return x.value


def to_class(c: Type[T], x: Any) -> dict:
    assert isinstance(x, c)
    return cast(Any, x).to_dict()


class Edge:
    name: Optional[str]
    src: Optional[str]
    dst: Optional[str]

    def __init__(self, name: Optional[str], src: Optional[str], dst: Optional[str]) -> None:
        self.name = name
        self.src = src
        self.dst = dst

    @staticmethod
    def from_dict(obj: Any) -> 'Edge':
        assert isinstance(obj, dict)
        name = from_union([from_str, from_none], obj.get("name"))
        src = from_union([from_str, from_none], obj.get("src"))
        dst = from_union([from_str, from_none], obj.get("dst"))
        return Edge(name, src, dst)

    def to_dict(self) -> dict:
        result: dict = {}
        result["name"] = from_union([from_str, from_none], self.name)
        result["src"] = from_union([from_str, from_none], self.src)
        result["dst"] = from_union([from_str, from_none], self.dst)
        return result


class Part(Enum):
    A1 = "A1"
    AMT_R = "AmtR"
    BYDV_J = "BydvJ"
    B_BA_B0064_RBS = "BBa_B0064_rbs"
    ECK120010876 = "ECK120010876"
    ECK120015440 = "ECK120015440"
    ECK120033736 = "ECK120033736"
    ECK120033737 = "ECK120033737"
    H1 = "H1"
    HLY_IIR = "HlyIIR"
    L3_S2_P21_TERMINATOR = "L3S2P21_terminator"
    L3_S2_P55 = "L3S2P55"
    LAC_I = "LacI"
    LMR_A = "LmrA"
    N1 = "N1"
    P3 = "P3"
    PHL_F = "PhlF"
    PSR_A = "PsrA"
    P_AMT_R = "pAmtR"
    P_CONST = "pCONST"
    P_HLY_IIR = "pHlyIIR"
    P_LMR_A = "pLmrA"
    P_PHL_F = "pPhlF"
    P_PSR_A = "pPsrA"
    P_TAC = "pTac"
    R1 = "R1"
    RIBO_J = "RiboJ"
    RIBO_J51 = "RiboJ51"
    RIBO_J53 = "RiboJ53"
    RIBO_J64 = "RiboJ64"
    SCM_J = "ScmJ"
    YFP = "YFP"


class Placement:
    position: Optional[int]
    direction: Optional[int]
    parts: Optional[List[Part]]

    def __init__(self, position: Optional[int], direction: Optional[int], parts: Optional[List[Part]]) -> None:
        self.position = position
        self.direction = direction
        self.parts = parts

    @staticmethod
    def from_dict(obj: Any) -> 'Placement':
        assert isinstance(obj, dict)
        position = from_union([from_int, from_none], obj.get("position"))
        direction = from_union([from_int, from_none], obj.get("direction"))
        parts = from_union([lambda x: from_list(Part, x), from_none], obj.get("parts"))
        return Placement(position, direction, parts)

    def to_dict(self) -> dict:
        result: dict = {}
        result["position"] = from_union([from_int, from_none], self.position)
        result["direction"] = from_union([from_int, from_none], self.direction)
        result["parts"] = from_union([lambda x: from_list(lambda x: to_enum(Part, x), x), from_none], self.parts)
        return result


class Node:
    name: Optional[str]
    node_type: Optional[str]
    partition_id: Optional[int]
    gate_type: Optional[str]
    placements: Optional[List[Placement]]

    def __init__(self, name: Optional[str], node_type: Optional[str], partition_id: Optional[int], gate_type: Optional[str], placements: Optional[List[Placement]]) -> None:
        self.name = name
        self.node_type = node_type
        self.partition_id = partition_id
        self.gate_type = gate_type
        self.placements = placements

    @staticmethod
    def from_dict(obj: Any) -> 'Node':
        assert isinstance(obj, dict)
        name = from_union([from_str, from_none], obj.get("name"))
        node_type = from_union([from_str, from_none], obj.get("nodeType"))
        partition_id = from_union([from_int, from_none], obj.get("partitionID"))
        gate_type = from_union([from_str, from_none], obj.get("gateType"))
        placements = from_union([lambda x: from_list(Placement.from_dict, x), from_none], obj.get("placements"))
        return Node(name, node_type, partition_id, gate_type, placements)

    def to_dict(self) -> dict:
        result: dict = {}
        result["name"] = from_union([from_str, from_none], self.name)
        result["nodeType"] = from_union([from_str, from_none], self.node_type)
        result["partitionID"] = from_union([from_int, from_none], self.partition_id)
        result["gateType"] = from_union([from_str, from_none], self.gate_type)
        result["placements"] = from_union([lambda x: from_list(lambda x: to_class(Placement, x), x), from_none], self.placements)
        return result


class Welcome:
    name: Optional[str]
    input_filename: Optional[str]
    nodes: Optional[List[Node]]
    edges: Optional[List[Edge]]

    def __init__(self, name: Optional[str], input_filename: Optional[str], nodes: Optional[List[Node]], edges: Optional[List[Edge]]) -> None:
        self.name = name
        self.input_filename = input_filename
        self.nodes = nodes
        self.edges = edges

    @staticmethod
    def from_dict(obj: Any) -> 'Welcome':
        assert isinstance(obj, dict)
        name = from_union([from_str, from_none], obj.get("name"))
        input_filename = from_union([from_str, from_none], obj.get("inputFilename"))
        nodes = from_union([lambda x: from_list(Node.from_dict, x), from_none], obj.get("nodes"))
        edges = from_union([lambda x: from_list(Edge.from_dict, x), from_none], obj.get("edges"))
        return Welcome(name, input_filename, nodes, edges)

    def to_dict(self) -> dict:
        result: dict = {}
        result["name"] = from_union([from_str, from_none], self.name)
        result["inputFilename"] = from_union([from_str, from_none], self.input_filename)
        result["nodes"] = from_union([lambda x: from_list(lambda x: to_class(Node, x), x), from_none], self.nodes)
        result["edges"] = from_union([lambda x: from_list(lambda x: to_class(Edge, x), x), from_none], self.edges)
        return result


def welcome_from_dict(s: Any) -> Welcome:
    return Welcome.from_dict(s)


def welcome_to_dict(x: Welcome) -> Any:
    return to_class(Welcome, x)


with open('/Users/chentianchi/Downloads/document.json') as f:
    data = json.load(f)

result = welcome_from_dict(data)

re_dict = result.to_dict()
