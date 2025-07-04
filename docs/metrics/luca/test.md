<!DOCTYPE html>
<html lang="fr">
<head>
  <meta charset="UTF-8">
  <title>Tableau Récapitulatif - Complet</title>
  <style>
    body {
      font-family: "Inter", sans-serif;
      margin: 20px;
      background-color: #f4f7f6;
      color: #333;
    }
    h1 {
      color: #2c3e50;
      text-align: center;
      margin-bottom: 30px;
    }
    table {
      width: 100%;
      border-collapse: collapse;
      margin: 0 auto;
      box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
      border-radius: 8px;
      overflow: hidden;
      table-layout: fixed;
    }
    th, td {
      padding: 10px 15px;
      border: 1px solid #ddd;
      vertical-align: top;
      line-height: 1.5;
      white-space: normal;
      word-wrap: break-word;
      word-break: break-word;
    }
    th {
      background-color: #4CAF50;
      color: white;
    }
    tr:nth-child(even) {
      background-color: #f2f2f2;
    }
    tr:hover {
      background-color: #e0e0e0;
    }
    td a {
      color: #007bff;
      text-decoration: none;
    }
    td a:hover {
      text-decoration: underline;
    }
    td span {
      display: inline-block;
      padding: 3px 8px;
      border-radius: 5px;
      font-size: 0.9em;
      font-weight: bold;
      text-align: center;
    }
    .batch-integration {
      background-color: #b10202;
      color: #ffcfc9;
    }
    .clustering {
      background-color: #e6cff2;
      color: #5a3286;
    }
    .denoising {
      background-color: #d4edbc;
      color: #11734b;
    }
    .dimension-reduction {
      background-color: #bfe1f6;
      color: #0a53a8;
    }
    .spatial-decomposition {
      background-color: #e8eaed;
      color: #434343;
    }
    .spatially-variable-genes {
      background-color: #c6dbe1;
      color: #215a6c;
    }
    .specificity {
      background-color: #ffe5a0;
      color: #473821;
    }
  </style>
</head>
<body>
  <h1>SUMMARY</h1>
  <table>
    <thead>
      <tr>
        <th>#</th>
        <th style="width:4cm;">Name</th>
        <th style="width:4cm;">Source</th>
        <th style="width:4cm;">Upload</th>
        <th style="width:3cm;">Date</th>
        <th style="width:4cm;">Name on GitHub</th>
        <th style="width:3cm;">Verification</th>
        <th style="width:3cm;">Check date</th>
        <th style="width:6cm;">Type</th>
        <th style="width:20cm;">Summary</th>
        <th style="width:6cm;">Interval</th>
      </tr>
    </thead>
    <tbody>
      <!-- Contenu ici : Copiez-collez ici les 44 lignes du tableau HTML que tu as déjà fourni -->
      <!-- Exemple de ligne (supprimer ce commentaire et coller le contenu complet) -->
      <tr>
        <td>1</td>
        <td>Adjusted Rand Index (ARI)</td>
        <td>Luecken et. al // OP</td>
        <td></td>
        <td>20/06/25</td>
        <td><a target="_blank" href="http://adjusted_rand_index.md/">adjusted_rand_index.md</a></td>
        <td></td>
        <td></td>
        <td><span class="batch-integration">Batch integration</span></td>
        <td>Adjusted Rand Index compares clustering overlap, correcting for random labels and considering correct overlaps and disagreements.</td>
        <td>[-1,1] the higher the better</td>
      </tr>
      <!-- etc. -->
    </tbody>
  </table>
</body>
</html>
