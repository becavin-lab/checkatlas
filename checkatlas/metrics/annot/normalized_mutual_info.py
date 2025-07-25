from sklearn.metrics import normalized_mutual_info_score


def run(annotation, ref_annotation):

    return normalized_mutual_info_score(annotation, ref_annotation)
