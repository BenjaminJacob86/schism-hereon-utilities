#!/usr/bin/env bash

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

. ""${SCRIPTPATH}"/read-secret.sh" --source-only

BASE_URL="https://gitlab.mercator-ocean.fr/api/v4/projects/2064/packages/helm/"

update_or_publish_chart() {
    CHART_NAME=$1
    echo "Pushing $CHART_NAME"

    if [[ $(helm plugin list | grep cm-push | wc -c) -eq 0 ]]; then
        helm plugin install https://github.com/chartmuseum/helm-push
    fi
    helm dependency update $CHART_NAME
    helm cm-push $CHART_NAME process-playground
}

get_only_modified_charts() {
    pushd "$SCRIPTPATH" > /dev/null
    LAST_MODIFIED_FILES=$(git diff origin/main...HEAD --name-only)
    for chart in *; do
        if [ -d "$chart" ]; then
            for modified_file in $LAST_MODIFIED_FILES; do
                if [[ $modified_file == ${chart}* ]]; then
                    echo $chart
                    break
                fi
            done
        fi
    done
    popd > /dev/null
}


publish_all_charts() {
    exportSecret "DEPLOY_TOKEN_USERNAME" "DEPLOY_TOKEN_USERNAME"
    exportSecret "DEPLOY_TOKEN_PASSWORD" "DEPLOY_TOKEN_PASSWORD"
    EXIT_CODE=0
    HELM_ENDPOINT="$BASE_URL/process-playground"
    helm repo add --username $DEPLOY_TOKEN_USERNAME --password $DEPLOY_TOKEN_PASSWORD process-playground $HELM_ENDPOINT
    pushd "$SCRIPTPATH" > /dev/null
    for chart in $(get_only_modified_charts); do
        update_or_publish_chart $chart || EXIT_CODE=1
    done
    popd > /dev/null
    exit $EXIT_CODE
}

lint_all_charts() {
    EXIT_CODE=0
    pushd "$SCRIPTPATH" > /dev/null
    for chart in $(get_only_modified_charts); do
        helm dependencies update $chart
        helm lint $chart --quiet || EXIT_CODE=1
    done
    popd > /dev/null
    exit $EXIT_CODE
}

$*
